classdef Mobility < StateFunction
    properties
    end
    
    methods
        function gp = Mobility(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'RelativePermeability', 'Viscosity'});
        end
        function mob = evaluateOnDomain(prop, model, state)
            [mu, kr] = prop.getEvaluatedDependencies(state, 'Viscosity', 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            if isfield(model.fluid, 'tranMultR')
                % Pressure dependent mobility multiplier 
                p = model.getProp(state, 'pressure');
                mult = model.fluid.tranMultR(p);
                mob = cellfun(@(x) x.*mult, mob, 'UniformOutput', false);
            end
            mv = cellfun(@(x) min(value(x)), mob);
            if model.verbose > 1 && any(mv < 0)
            	warning('Negative mobilities detected! Capping to zero.')
                mob = cellfun(@(x) max(x, 0), mob, 'UniformOutput', false);
            end
        end
    end
end