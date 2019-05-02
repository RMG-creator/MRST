classdef SurfactantFlowPropertyFunctions < BlackOilFlowPropertyFunctions
    properties
        CapillaryNumber
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@BlackOilFlowPropertyFunctions(model);
            props.RelativePermeability = SurfactantRelativePermeability(model);
            props.CapillaryNumber = CapillaryNumber(model);
        end
    end
end