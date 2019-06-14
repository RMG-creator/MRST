classdef BlackOilSurfactantViscosity < BlackOilViscosity
    methods
        function gp = BlackOilSurfactantViscosity(prop, varargin)
            gp@BlackOilViscosity(prop, varargin{:});
            gp = addPropertyDependence(gp, {'SurfactantViscMultiplier'});
            gp = gp.dependsOn('surfactant', 'state');
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            cs = model.getProps(state, 'surfactant');
            muWMults = prop.getEvaluatedDependencies(state, 'SurfactantViscMultiplier');
            mu = prop.evaluateOnDomain@BlackOilViscosity(model, state);
            mu{1} = model.fluid.muWSft(cs);
            mu{1} = mu{1}.*muWMults;            
        end
    end
end