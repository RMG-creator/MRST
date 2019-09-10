classdef PhasePressures < StateFunction
    properties
    end
    
    methods
        function gp = PhasePressures(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'CapillaryPressure', 'Pressure'});
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            [p, pc] = prop.getEvaluatedDependencies(state, 'Pressure', 'CapillaryPressure');
            nph = numel(pc);
            p_phase = cell(1, nph);
            for i = 1:nph
                if isempty(pc{i})
                    p_phase{i} = p;
                else
                    p_phase{i} = p + pc{i};
                end
            end
        end
    end
end