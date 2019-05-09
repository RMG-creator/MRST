classdef PolymerFlowPropertyFunctions < BlackOilFlowPropertyFunctions

    properties
        PolymerAdsorption;
        PolymerViscMultiplier;
        EffectiveMixturePolymerViscMultiplier;
    end
    
    methods
        function props = PolymerFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.PolymerViscMultiplier = PolymerViscMultiplier(model, sat);
            props.EffectiveMixturePolymerViscMultiplier = EffectiveMixturePolymerViscMultiplier(model, sat);
            props.Viscosity             = PolymerViscosity(model, sat);
            props.PolymerAdsorption     = PolymerAdsorption(model, sat);
        end
    end
end