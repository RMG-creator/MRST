classdef PressureNaturalVariablesModel < NaturalVariablesCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        useIncTolPressure
        singlePhaseStrategy = 'numerical';
        twoPhaseStrategy = 'numerical';
        singlePhaseDifferentiation = 'numerical';
        twoPhaseDifferentiation = 'numerical';
    end
    
    methods
        function model = PressureNaturalVariablesModel(G, rock, fluid, compFluid, varargin)
            
            model = model@NaturalVariablesCompositionalModel(G, rock, fluid, compFluid);
            model.useIncTolPressure = true;
            model.useIncTolComposition = true;
            model.allowLargeSaturations = true;
            model.maxPhaseChangesNonLinear = 20;
            model = merge_options(model, varargin{:});
            model.outputFluxes = true;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsNaturalVariables(state0, state, model, dt, ...
                            drivingForces, 'reduceToPressure', true, varargin{:});
            problem = problem.assembleSystem();
            
            state.pRes = norm(problem.b, inf);
            state.b = problem.b;
            
            
            state.w = problem.w;
            state.w_p = state.pressure;
            if problem.iterationNo == 1
                state.z0 = state.components;
                if model.water
                   mass = state.s.*state.rho;
                   state.sM = mass(:, 1)./sum(mass, 2);
                   state.sW0 = state.s(:, 1);
                end
                
                problem.state = state;
            end
        end
        function [state, report] = updateState(model, state, problem, dv, drivingForces)
            
            state.dpressure = dv{1};
            [state, report] = updateState@NaturalVariablesCompositionalModel(model, state, problem, dv, drivingForces);
        end
        
        function  [values, tolerances, names] = getConvergenceValues(model, problem)
            [values, tolerances, names] = getConvergenceValues@NaturalVariablesCompositionalModel(model, problem);
            if model.water
                isWat = strcmpi(names, 'water');
                tolerances = tolerances(~isWat);
                names = names(~isWat);
                values = values(~isWat);
            end
            if ~model.useIncTolPressure
                sub = strcmpi(names, 'deltap');
                values(sub) = norm(problem.b(1:model.G.cells.num), inf);
                tolerances(sub) = model.nonlinearTolerance;
                names{1} = 'Pressure';
            end
        end
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@NaturalVariablesCompositionalModel(model, state0, state, dt, drivingForces);
            state.switchCount_p = state.switchCount;
        end
    end
end