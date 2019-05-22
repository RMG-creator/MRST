classdef PhasePressures < GridProperty
    properties
    end
    
    methods
        function gp = PhasePressures(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'CapillaryPressure'});
            gp = gp.dependsOn({'pressure'}, 'state');
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
%             fluid = model.fluid;
            p = model.getProps(state, 'Pressure');
            pc = prop.getEvaluatedDependencies(state, 'CapillaryPressure');
%             if model.surfactant
%                 pc{1} = pc{1}.*fluid.ift(cs)/fluid.ift(0);
%                 pc{2} = pc{2}.*fluid.ift(cs)/fluid.ift(0);
%             end
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