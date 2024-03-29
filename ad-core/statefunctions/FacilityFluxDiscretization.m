classdef FacilityFluxDiscretization < StateFunctionGrouping
    properties
        PhaseFlux
        ComponentTotalFlux
        ComponentPhaseFlux
        PerforationPressureGradient
        WellIndex
        FacilityWellMapping
    end
    
    methods
        function props = FacilityFluxDiscretization(model)
            props.structName = 'FacilityFluxProps';
            
            props.PhaseFlux = WellPhaseFlux(model);
            props.ComponentTotalFlux = ComponentTotalFlux(model);
            props.ComponentPhaseFlux = WellComponentPhaseFlux(model);
            props.PerforationPressureGradient = PerforationPressureGradient(model);
            props.WellIndex = WellIndex(model);
            props.FacilityWellMapping = FacilityWellMapping(model);
        end
    end
end