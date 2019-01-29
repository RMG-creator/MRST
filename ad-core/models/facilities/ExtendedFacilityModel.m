classdef ExtendedFacilityModel < FacilityModel
    properties
        
    end
    
    methods
        
        function n = getNumberOfComponents(fm)
            n = fm.ReservoirModel.getNumberOfComponents();
        end
        
        function n = getNumberOfPhases(fm)
            n = fm.ReservoirModel.getNumberOfPhases();
        end
        
        function src = getComponentSourceTerm(fm, state, name)
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@FacilityModel(model, state, state0, dt, drivingForces);
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);
            
            if nw > 0
                rho = model.ReservoirModel.getProps(state, 'Density');
                rho = [rho{:}];
                [wc, perf2well] = model.getActiveWellCells(wellSol);
                rho = rho(wc, :);
                for i = 1:nw
                    wellNo = actWellIx(i);
                    wm = model.WellModels{wellNo};
                    % Possible index error for i here - should it be
                    % wellno?
                    rho_i = rho(perf2well == i, :);
                    wellSol(wellNo) = wm.updateConnectionPressureDropState(model.ReservoirModel, wellSol(wellNo), rho_i, rho_i);
                end
            end
        end
        
        function containers = getPropertyFunctions(model)
            containers = getPropertyFunctions@PhysicalModel(model);
            assert(not(isempty(model.FacilityFluxDiscretization)), ...
                'PropertyFunctions not initialized - did you call "validateModel"?');
            containers = [containers, {model.FacilityFluxDiscretization}];
        end
    end
end