classdef MEORXModel < TwoPhaseOilWaterModel
    % Oil/water/microbial system
    % This model is a two phase oil/water model, extended with the 
    % microbial phase in addition.
    % Microbe effects currently available for simulation:
    % Nothing
    
    properties
        % Microbes and Nutrients present
        microbe_0
        microbe_1
        nutrient
        pollutes
        biofilm
        yield
        growth
        biosurf
        biopoly
    end
    
    methods
        function model = MEORXModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G,rock, fluid);
            
            model.microbe_0 = true;
            model.microbe_1 = true;
            model.pollutes = true;
            model.nutrient = true;
            model.yield = true;
            model.growth = true;
            % type/types of metabolite to be defined manually
            model.biosurf = false;
            model.biopoly = false;
            % not available yet
            model.outputFluxes = false;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWMEOR', 'bhp'};
            
            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state,...
                dt, drivingForces, varargin)
            [problem, state] = equationsMEORX(state0, state, model,...
                dt, drivingForces, varargin{:});
        end
        
        function [fn, index] = getVariableField(model, name)
            % add metabolites here
            switch(lower(name))
                case 'microbe_0'
                    fn = 'm_0';
                    index = 1;
                case 'microbe_1'
                    fn = 'm_1';
                    index = 1;
                case 'nutrient'
                    fn = 'n';
                    index = 1;
                case 'pollutes'
                    fn = 'pol';
                    index = 1;
                case 'biofilm'
                    fn = 'bio';
                    index = 1;
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(model, name);
            end
        end
    end
end