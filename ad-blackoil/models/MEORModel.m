classdef MEORModel < TwoPhaseOilWaterModel
    % Oil/water/microbial system
    % This model is a two phase oil/water model, extended with the 
    % microbial phase in addition.
    % Microbe effects currently available for simulation:
    % Nothing
    
    properties
        % Microbes and Nutrients present
        microbe
        nutrient
    end
    
    methods
        function model = MEORModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G,rock, fluid);
            
            model.microbe = true;
            model.nutrient = true;
            % not available yet
            model.outputFluxes = false;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWMEOR', 'bhp'};
            
            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state,...
                dt, drivingForces, varargin)
            [problem, state] = equationsMEOR(state0, state, model,...
                dt, drivingForces, varargin{:});
        end
        
        function [fn, index] = getVariableField(model, name)
            % add metabolites here
            switch(lower(name))
                case 'microbe'
                    fn = 'm';
                    index = 1;
                case 'nutrient'
                    fn = 'n';
                    index = 1;
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(model, name);
            end
        end
    end
end