classdef PolymerFlowPropertyFunctions < BlackOilFlowPropertyFunctions

    properties
        PolymerAdsorption;
        PolymerViscMultiplier;
        PolymerPhaseFlux;
        FaceConcentration;
    end
    
    methods
        function props = PolymerFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.PolymerViscMultiplier = PolymerViscMultiplier(model, sat);
            props.Viscosity             = PolymerViscosity(model, sat);
            props.PolymerAdsorption     = PolymerAdsorption(model, sat);
            props.FaceConcentration     = FaceConcentration(model, sat);
            props.PolymerPhaseFlux      = PolymerPhaseFlux(model, sat);
        end
    end
end