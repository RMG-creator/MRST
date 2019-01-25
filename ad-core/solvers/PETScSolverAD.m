classdef PETScSolverAD < LinearSolverAD
    % Linear solver that calls a PETSc solver
    %
    % SYNOPSIS:
    %   solver = PETScSolverAD()
    %
    % DESCRIPTION:
    %    AD-interface for the PETSc interface.
    %
    % NOTE:
    %    This solver requires PETSc to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`
    
    properties
        petsc_opt
    end
    
    methods
        function solver = PETScSolverAD(varargin)
            require linearsolvers
            solver = solver@LinearSolverAD();
            
            % Default options: FIXME check all other options (tolerance,
            % maxIterations, extraReport etc)
            
            solver = merge_options(solver, varargin{:});
        end
        
        function [result, report] = solveLinearSystem(solver, A, b)
            solver = merge_petsc_opt(solver);
            
            [result, stats] = callPETSc(A, b, 'petsc_opt', solver.petsc_opt);
            
            report = struct('converged',  stats.reason > 0, ...
                'residual',   stats.err,...
                'iterations', stats.nIter);
            
            % Debugging
            report
            x = A\b;
            plot(result-x), title 'petsc - mldivide', drawnow
        end 
        
        function solver = merge_petsc_opt(solver)
            % Append options to petsc. Note that solver.tolerance and
            % solver.maxIterations are set for the ksp solver which may or
            % may not be used by petsc. Specific petsc solver parameters
            % should be set in the petsc_opt string.
            solver.petsc_opt = [solver.petsc_opt, ' -ksp_rtol ', num2str(solver.tolerance)];
            solver.petsc_opt = [solver.petsc_opt, ' -ksp_max_it ', num2str(solver.maxIterations)];
        end
    end
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
