classdef SurfactantFlowPropertyFunctions < BlackOilFlowPropertyFunctions
    properties
        CapillaryNumber;
        SurfactantViscMultiplier;
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.RelativePermeability = SurfactantRelativePermeability(model);
            props.CapillaryNumber = CapillaryNumber(model);
            props.SurfactantViscMultiplier = SurfactantViscMultiplier(model, sat);
            props.Viscosity = BlackOilSurfactantViscosity(model, sat);
            props.PhasePressures = SurfactantPhasePressures(model, sat);
        end
    end
end