classdef BlackOilPolymerViscosity < BlackOilViscosity
    methods
        function gp = BlackOilPolymerViscosity(prop, varargin)
            gp@BlackOilViscosity(prop, varargin{:});
            gp = addPropertyDependence(gp, {'PolymerViscMultiplier'});
            gp = addPropertyDependence(gp, {'polymer'}, 'state');
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            muWMult = prop.getEvaluatedDependencies(state, 'PolymerViscMultiplier');
            mu = prop.evaluateOnDomain@BlackOilViscosity(model, state);
            mu{1} = mu{1}.*muWMult;            
        end
    end
end