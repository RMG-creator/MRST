classdef BlackOilFluidModel < ThreePhaseBlackOilModel

    properties
        primaryVarNames;
    end

    methods
        function model = BlackOilFluidModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);

            model.disgas          = true;
            model.vapoil          = false;
            model.primaryVarNames = {'pressure', 'sW', 'x'}; % well variables
                                                             % not included
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The BlackOilFluidModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function varnames = getAllPrimaryVariables(model, state)
            varnames = model.primaryVarNames;
            [~, wellVarNames, ~] = ...
                model.FacilityModel.getAllPrimaryVariables(state.wellSol);
            varnames = {varnames{:}, wellVarNames{:}};
        end

        function fds = getAllVarsNames(model)
            fds = {'wellSol', 'pressure', 's', 'rs', 'rv', 'sW', 'sG', 'sO', ...
                   'water', 'oil', 'gas'};
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            vars = problem.primaryVariables;
            ind = false(size(vars));
            fluidvars = model.getAllPrimaryVariables(state);

            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, ...
                                                              problem, dx, ...
                                                              drivingForces);

        end


    end

end