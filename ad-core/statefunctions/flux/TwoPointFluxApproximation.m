classdef TwoPointFluxApproximation < PermeabilityGradientDiscretization
    properties
    end
    
    methods
        function tpfa = TwoPointFluxApproximation(model)
            assert(isfield(model.operators, 'T'), 'Transmissibility must be present in model.operators for TPFA discretization.');
        end

        function v = getPermeabilityGradient(tpfa, model, state, dp)
            T = model.getProp(state, 'transmissibility');
            v = T.*dp;
        end
    end
end