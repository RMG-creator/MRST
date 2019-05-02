classdef PolymerFlowPropertyFunctions < BlackOilFlowPropertyFunctions

    properties
        BlackoilMobility;
        PolymerAdsorption;
        PolymerViscMultiplier;
    end
    
    methods
        function props = PolymerFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.Mobility              = PolymerMobility(model, sat);
            props.PolymerViscMultiplier = PolymerViscMultiplier(model, sat);
            props.BlackoilMobility      = Mobility(model, sat);
            props.PolymerAdsorption     = PolymerAdsorption(model, sat);
        end
    end
end