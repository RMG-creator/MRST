classdef ThreePhaseCompositionalModel < ReservoirModel
% Compositional model with up to two phases + optional aqua phase
    properties
        EOSModel
        EOSNonLinearSolver
        dzMaxAbs
        incTolPressure
    end
    
    methods
        function model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model.EOSModel = EquationOfStateModel(G, compFluid);
            
            model.EOSModel.verbose = false;
            
            model.EOSNonLinearSolver = NonLinearSolver();
            model.EOSNonLinearSolver.verbose = false;
            model.EOSNonLinearSolver.maxIterations = 100;
            model.EOSNonLinearSolver.maxTimestepCuts = 0;
            model.EOSNonLinearSolver.useRelaxation = true;
            model.EOSNonLinearSolver.continueOnFailure = true;
            model.EOSNonLinearSolver.errorOnFailure = false;
            
            model.nonlinearTolerance = 0.01;
            model.incTolPressure = 1e-3;
            
            model.water = true;
            model.oil = true;
            model.gas = true;
            
            model.dzMaxAbs = 0.1;
            model.dsMaxAbs = 0.1;
            
            model.minimumPressure = 0;
            model = merge_options(model, varargin{:});
            
            if model.water
                model.saturationVarNames = {'sw', 'so', 'sg'};
                model.wellVarNames = {'qWs', 'qOs', 'qGs', 'bhp'};
            else
                model.saturationVarNames = {'so', 'sg'};
                model.wellVarNames = {'qOs', 'qGs', 'bhp'};
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsCompositional(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
        
        function [fn, index] = getVariableField(model, name)
            switch(lower(name))
                case {'z', 'components'}
                    % Overall mole fraction
                    fn = 'components';
                    index = ':';
                case {'x', 'liquidmf'}
                    % Liquid phase mole fraction
                    fn = 'x';
                    index = ':';
                case {'y', 'vapormf'}
                    % Vapor phase mole fraction
                    fn = 'y';
                    index = ':';
                case {'l', 'liquid'}
                    % Liquid mole fraction
                    fn = 'L';
                    index = 1;
                case {'z_l'}
                    % Liquid compressibility
                    fn = 'Z_L';
                    index = 1;
                case {'z_v'}
                    % Vapor compressibility
                    fn = 'Z_V';
                    index = 1;
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state0 = state;
            
            var0 = problem.primaryVariables;
            vars = var0;
            removed = false(size(vars));

            % Saturation update
            wIx = strcmpi(vars, 'sw');
            
            if any(wIx)
                state = model.updateStateFromIncrement(state, dx{wIx}, problem, 'sW', model.dsMaxRel, model.dsMaxAbs);
                state.s(:, 1) = max(min(state.s(:, 1), 1), 0);
                
                state.s(:, 2) = state.Z_L.*state.L./(state.Z_L.*state.L + state.Z_V.*(1 - state.L));
                state.s(:, 3) = 1 - state.s(:, 2) - state.s(:, 3);

                [vars, ix] = model.stripVars(vars, {'sw'});
                removed(~removed) = removed(~removed) | ix;
            end
            
            % Components
            cnames = model.EOSModel.fluid.names;
            ncomp = numel(cnames);
            ok = false(ncomp, 1);

            z = state.components;
            rm = 0;
            
            if 0
                dz = cell(ncomp, 1);
                
                rm = 0;
                for i = 1:ncomp
                    name = cnames{i};
                    cix = strcmpi(var0, name);
                    if ~any(cix)
                        continue
                    end
                    
                    dz{i}= dx{cix};
                    ok(i) = true;
                    rm = rm - dz{i};
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;

                end
                
                if any(ok)
                    dz{~ok} = rm;
                    w = ones(model.G.cells.num, 1);
                    for i = 1:ncomp
                        dzcap = (min(max(z{i} + dz{i}, 0), 1) - z{i});
                        w = min(min(dzcap./dz{i}, 1), w);
                        w(~isfinite(w)) = 0;
                    end

                    
                    for i = 1:ncomp
                        z{i} = z{i} + w.*dz{i};
                    end
                end
                if any(ok)
                    % We had components as active variables somehow
                    assert(nnz(~ok) == 1)
                    z{~ok} = min(max(z{~ok} + rm, 0), 1);
                    sz = sum([z{:}], 2);
                    z = cellfun(@(x) x./sz, z, 'UniformOutput', false);
                    state.components = z;

                    v  = 1 - state.s(:, 1);
                    v0 = 1 - state0.s(:, 1);
                    state.dz = computeChange(state.components, state0.components, v, v0);
                end
            else
                for i = 1:ncomp
                    name = lower(cnames{i});
                    cix = strcmpi(var0, name);
                    if any(cix)
                        z0 = z{i};

                        dz = dx{cix};
                        if isfinite(model.dzMaxAbs)
                            dz = sign(dz).*min(abs(dz), model.dzMaxAbs);
                        end
                        z{i} = min(max(z0 + dz, 0), 1);

                        ok(i) = true;
                        [vars, ix] = model.stripVars(vars, {name});
                        removed(~removed) = removed(~removed) | ix;

                        rm = rm - (z{i} - z0);
                    end
                end
                if any(ok)
                    % We had components as active variables somehow
                    assert(nnz(~ok) == 1)
                    z{~ok} = min(max(z{~ok} + rm, 0), 1);
                    sz = sum([z{:}], 2);
                    z = cellfun(@(x) x./sz, z, 'UniformOutput', false);
                    state.components = z;
                    if model.water
                        v  = 1 - state.s(:, 1);
                        v0 = 1 - state0.s(:, 1);
                    else
                        [v, v0] = deal(1);
                    end
                    state.dz = computeChange(state.components, state0.components, v, v0);
                end
                
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dx(removed) = [];
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);

            p0 = state0.pressure;
            range = max(p0) - min(p0);
            if range == 0
                range = 1;
            end
            state.dpRel = (state.pressure - p0)./range;
            state.dpAbs = state.pressure - p0;
            
            
        end
        
        function state = computeFlash(model, state, dt, iteration)
            if nargin < 4
                iteration = 1;
            end
            transportProblem = isa(model, 'TransportCompositionalModel');
            if iteration == 1 && ~transportProblem
                state.eos.iterations = 0;
                state.eos.cellcount = 0;
            end
            
            if isfield(state, 'L');
                L0 = state.L;
            else
                L0 = [];
            end
            [state, report] = model.EOSNonLinearSolver.solveTimestep(state, dt, model.EOSModel);
            
            % stability chop for liquid fraction

            if ~isempty(report)
%                 fprintf('Flash computed in %d iterations\n', report.StepReports{1}.Iterations);
                state.eos.iterations = state.eos.iterations + report.StepReports{1}.Iterations;
                state.eos.cellcount = state.eos.cellcount + sum(cellfun(@(x) x.ActiveCells, report.StepReports{1}.NonlinearReport));
                
                if ~report.StepReports{1}.Converged
                    disp('!!! Flash did not converge');
                    fprintf('Final residuals: ')
                    fprintf('%1.4g ', report.StepReports{1}.NonlinearReport{end}.Residuals)
                    fprintf('after %d iterations\n', report.StepReports{1}.Iterations);
                end
            end
            L = state.L;
            Z_L = state.Z_L;
            Z_V = state.Z_V;
            
            if model.water
                sW = state.s(:, 1);
                sL = (1 - sW).*L.*Z_L./(L.*Z_L + (1-L).*Z_V);
                sV = 1 - sW - sL;
                state.s = [sW, sL, sV];
            else
                sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
                sV = 1 - sL;
                state.s = [sL, sV];
            end
            if ~isempty(L0) && 0
                tol = 1e-4;
                state.L(L == 1 & L0 < 1 - tol) = 1 - tol;
                state.L(L0 == 1 & L < 1 - tol) = 1 - tol;
                state.L(L == 0 & L0 > tol) = tol;
                state.L(L > tol & L0 == 0) = tol;
            end
        end
    
        function [xM, yM, sL, sV, rhoL, rhoV, muL, muV] = computeTwoPhaseFlowProps(model, state, p, temp, z)
            isADI = isa(p, 'ADI') || isa(temp, 'ADI') || any(cellfun(@(x) isa(x, 'ADI'), z));
            
            if isADI
                [x, y, L] = model.EOSModel.getPhaseFractionAsADI(state, p, temp, z);
                [Z_L, Z_V] = model.EOSModel.getCompressibility(state, p, temp, x, y);
            else
                [x, y, L, Z_L, Z_V] = model.getProps(state, 'x', 'y', 'L', 'Z_L', 'Z_V');
            end
            xM = model.EOSModel.getMassFraction(x);
            yM = model.EOSModel.getMassFraction(y);

            rhoL = model.EOSModel.computeDensity(p, x, Z_L, temp, true);
            rhoV = model.EOSModel.computeDensity(p, y, Z_V, temp, false);

            sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
            sV = 1 - sL;
            if nargout > 6
                muL = model.EOSModel.computeViscosity(p, rhoL, temp, x, true);
                muV = model.EOSModel.computeViscosity(p, rhoV, temp, y, false);
            end
        end

        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            if ~isa(model, 'PressureCompositionalModel')
                names = problem.equationNames;

                offset = model.water;
                ncomp = numel(model.EOSModel.fluid.names);
                if isa(model, 'TransportCompositionalModel')
                    % We might be missing a 
                    ncomp = ncomp - model.conserveWater;
                    offset = model.conserveWater;
                end

                values = cellfun(@(x) norm(double(x), inf), problem.equations);

                if problem.iterationNo > 1
                    v = problem.state.dz(1:ncomp);
                else
                    v = inf(1, ncomp);
                end

                values((1+offset):(ncomp+offset)) = v;

                [conv_wells, ~, isWell] = checkWellConvergence(model, problem);
                nonWellValues = values(~isWell);


                convergence = all(nonWellValues <= model.nonlinearTolerance) && all(conv_wells);
            else
                [convergence, values, names] = checkConvergence@ReservoirModel(model, problem, varargin{:});
            end 
            if ~isa(model, 'TransportCompositionalModel') && ~isa(model, 'PressureCompositionalModel');
                if isfield(problem.state, 'dpRel')
                    dp = norm(problem.state.dpRel, inf);
                else
                    dp = inf;
                end
                values = [dp, values];
                convergence = convergence && dp <= model.incTolPressure;
                names = ['deltaP', names];
            end
        end
    end
end

function dz = computeChange(z, z0, s, s0)
    z_prev = cellfun(@(x) x.*s0, z0, 'UniformOutput', false);
    z_curr = cellfun(@(x) x.*s, z, 'UniformOutput', false);

    dz = cellfun(@(x, y) norm(x - y, inf), z_curr, z_prev);
end