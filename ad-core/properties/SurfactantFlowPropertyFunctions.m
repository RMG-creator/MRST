classdef SurfactantFlowPropertyFunctions < BlackOilFlowPropertyFunctions
    properties
        CapillaryNumber
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            props.RelativePermeability = SurfRelativePermeability(model);
            props.CapillaryNumber = CapillaryNumber(model);
        end
    end
end