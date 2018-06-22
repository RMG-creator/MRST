classdef AMGCL_CPRSolverAD < AMGCLSolverAD
    % Linear solver that calls external compiled multigrid solver
    %
    % SYNOPSIS:
    %   solver = AMGCLSolverAD()
    %
    % DESCRIPTION:
    %    AD-interface for the AMGCL interface.
    %
    % NOTE:
    %    This solver requires AMGCL to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`

   properties
       doApplyScalingCPR
       trueIMPES % Use true impes decoupling strategy (if supported by model)
       useSYMRCMOrdering
   end
   methods
       function solver = AMGCL_CPRSolverAD(varargin)
            require linearsolvers
            solver = solver@AMGCLSolverAD();
            solver.trueIMPES    = false;
            solver.doApplyScalingCPR = true;
            solver.reduceToCell = true;
            solver.tolerance    = 1e-6;
            solver.useSYMRCMOrdering = true;
            
            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
           [result, report] = solver.callAMGCL_MEX(A, b, 2);
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           problem = solver.prepareProblemCPR(problem, model);
           [dx, result, report] = solveLinearProblem@LinearSolverAD(solver, problem, model);
       end
       
       function [grad, result, report] = solveAdjointProblem(solver, problemPrev,...
                problemCurr, adjVec, objective, model) %#ok
            % Solve an adjoint problem.
            timer = tic();
            problemCurr = problemCurr.assembleSystem();
            
            % Maybe formalize the control variables a bit in the future
            % sometime...
            if iscell(objective)
                objective = objective{:};
            end
            objective = combineEquations(objective);
            assert(isa(objective, 'ADI'), 'Objective function was not of type ADI.');
            b = -(objective.jac{1})';
            if ~isempty(adjVec)
                problemPrev = problemPrev.assembleSystem();
                b = b - problemPrev.A'*adjVec;
            end
            A = problemCurr.A;
            b = full(b);
            % Apply scaling
            %[A, b, scaling] = solver.applyScaling(A, b);
            % Reduce system (if requested)
            [A, b, lsys] = solver.reduceLinearSystemAdjoint(A, b);
            % Reorder linear system
            [A, b] = solver.reorderLinearSystem(A, b);
            %% make pressure like system
            scale = model.getScalingFactorsCPR(problemCurr, problemCurr.equationNames, solver);
            nc=size(A,1)/2  
            bindex=repmat(2,nc,1);
            %bb=[scale{1} scale{2};0 1];
            bb=[scale{1} scale{2};-scale{2} scale{1}];            
            [ind1, ind2] = blockDiagIndex(bindex, bindex);
            val=repmat(reshape(bb,[],1),nc,1);
            tt = sparse(ind1,ind2,val,size(A,1),size(A,2));
            A=tt*A;%b=tt*b;            
            
            if(false)
                sz = size(A);
                assert(sz(1) == sz(2), 'Matrix must be square!');
                n = sz(1);
                d = repmat([-barsa,1]',nc,1);
                %d=1./d;
                I = (1:n)';
                M = sparse(I, I, d, n, n);
            else
               M = solver.getDiagonalInverse(A); 
            end
            %[A, b, scaling] = solver.applyScaling(A, b);
            A=A*M;b=M*b;
            % Apply transpose
            A = A';
            t_prepare = toc(timer);
            % Solve system
            disp(max(b))
            [result, report] = solver.solveLinearSystem(A, b);            
            result=tt'*result;
            any(isnan(result))
            
            
            t_solve = toc(timer) - t_prepare;
            % Permute system back
            result = solver.deorderLinearSystemAdjoint(result);
            % Recover eliminated variables on linear level
            result = solver.recoverLinearSystemAdjoint(result, lsys);
            % Undo scaling
            %result = solver.undoScalingAdjoint(result, scaling);
            report.SolverTime = toc(timer);
            report.LinearSolutionTime = t_solve;
            report.preparationTime = t_prepare;
            report.postprocessTime = report.SolverTime - t_solve - t_prepare;
            grad = solver.storeIncrements(problemCurr, result);
       end
        

       function setSRelaxation(solver, v)
           solver.amgcl_setup.s_relaxation = translateOptionsAMGCL('relaxation', v);
       end
       
       function problem = prepareProblemCPR(solver, problem, model)
           n = model.G.cells.num;
           if solver.amgcl_setup.block_size == 0
               % Solver has not been told about block size, try to compute
               % it from what we are given.
               s = getSampleAD(problem.equations{:});
               nv = s.getNumVars();
               isCell = problem.indexOfType('cell');
               solver.amgcl_setup.block_size = sum(nv(isCell)/n);
           end
           
           % Get and apply scaling
           if solver.doApplyScalingCPR
               scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver);
               
               if solver.amgcl_setup.use_drs
                   % Solver will take the sum for us, we just weight each
                   % equation. Note: This is not the entirely correct way
                   % of doing this, as solver could do this by itself.
                   for i = 1:numel(scale)
                       if ~strcmpi(problem.types{i}, 'cell')
                           continue
                       end
                       if (numel(scale{i}) > 1 || scale{i} ~= 0)
                           problem.equations{i} = problem.equations{i}.*scale{i};
                       end
                   end
               else
                   % We form pressure equation for the solver.
                   e = 0;
                   for i = 1:numel(scale)
                       if ~strcmpi(problem.types{i}, 'cell')
                           continue
                       end
                       if (numel(scale{i}) > 1 || scale{i} ~= 0)
                           e = e + problem.equations{i}.*scale{i};
                       end
                   end
                   problem.equations{1} = e;
               end
           end
           m = solver.amgcl_setup.block_size;
           assert(m > 0);
           
           if isempty(solver.keepNumber)
               if solver.reduceToCell
                   % Will be reduced to ncell by block_size syste,
                   ndof = n*m;
               else
                   % We have no idea and should check
                   problem = problem.assembleSystem();
                   ndof = size(problem.A, 1);
                   if solver.amgcl_setup.active_rows == 0
                       % Only the first n*m entries are cell-wise
                       % variables, tell the solver this
                       solver.amgcl_setup.active_rows = n*m;
                   end
               end
           else
               ndof = solver.keepNumber;
           end
           
           
           if isempty(solver.variableOrdering) || numel(solver.variableOrdering) ~= ndof
               if solver.useSYMRCMOrdering
                   sym_ordering = getGridSYMRCMOrdering(model.G);
               else
                   sym_ordering = [];
               end
               ordering = getCellMajorReordering(n, m, 'ndof', ndof, 'cell_ordering', sym_ordering);
               solver.variableOrdering = ordering;
               if isempty(solver.equationOrdering) || numel(solver.equationOrdering) ~= ndof
                   solver.equationOrdering = ordering;
               end
           end
       end
   end
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
