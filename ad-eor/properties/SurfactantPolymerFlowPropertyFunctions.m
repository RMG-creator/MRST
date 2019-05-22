classdef SurfactantPolymerFlowPropertyFunctions < SurfactantFlowPropertyFunctions

    properties
        PolymerAdsorption;
        SurfactantAdsorption;
        PolymerViscMultiplier;
        SurfactantViscMultiplier;
        EffectiveMixturePolymerViscMultiplier;
        SurfactantPhasePressures;
    end
    
    methods
        function props = SurfactantPolymerFlowPropertyFunctions(model)
            props = props@SurfactantFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.PolymerViscMultiplier = PolymerViscMultiplier(model, sat);
            props.SurfactantViscMultiplier = SurfactantViscMultiplier(model, sat);
            props.EffectiveMixturePolymerViscMultiplier = EffectiveMixturePolymerViscMultiplier(model, sat);
            props.Viscosity             = SurfactantPolymerViscosity(model, sat);
            props.PolymerAdsorption     = PolymerAdsorption(model, sat);
            props.SurfactantAdsorption     = SurfactantAdsorption(model, sat);
            props.PhasePressures     = SurfactantPhasePressures(model, sat);
        end
    end
end