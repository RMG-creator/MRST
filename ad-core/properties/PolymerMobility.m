classdef PolymerMobility < GridProperty
    properties
    end

    methods
        function gp = PolymerMobility(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'BlackoilMobility', 'PolymerViscMultiplier'});
        end

        function mob = evaluateOnDomain(prop, model, state)
            [mob, muWMult] = prop.getEvaluatedDependencies(state, 'BlackoilMobility', 'PolymerViscMultiplier');
            mob{1} = mob{1}./muWMult; % get water index in the proper way
        end
    end
end