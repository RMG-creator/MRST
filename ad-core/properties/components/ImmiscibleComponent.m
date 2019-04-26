classdef ImmiscibleComponent < ComponentImplementation
    properties
        phaseIndex % Index of phase this component belongs to
    end
    
    methods
        function c = ImmiscibleComponent(name, phase)
            c@ComponentImplementation(name);
            c.phaseIndex = phase;
            c = c.dependsOn('Density');
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            c = getComponentDensity@ComponentImplementation(component, model, state, varargin{:});
            rho = model.getProp(state, 'Density');
            c{component.phaseIndex} = rho{component.phaseIndex};
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{component.phaseIndex} = 1;
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            c = component.getPhaseComposition(model, state);
        end
        
        function c = getPhaseComponentFractionWell(component, model, state, W)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(W.compi);
            index = component.phaseIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{index} = ci;
            end
        end
    end
end