classdef PolymerFluxDiscretization < FluxDiscretization
    properties
        PolymerPhaseFlux % Polymer phase volumetric fluxes
        FaceConcentration % Polymer or surfactant concentration on face
    end

    methods
        function props = PolymerFluxDiscretization(model)
            props = props@FluxDiscretization(model);
            if isempty(model.operators)
                N = getNeighbourship(model.G);
                nc = model.G.cells.num;
                nf = size(N, 1);
                up = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);
            else
                up = model.operators.faceUpstr;
            end
            upstr = UpstreamFunctionWrapper(up);
            props.PolymerPhaseFlux = PolymerPhaseFlux(model);            
            props.FaceConcentration = FaceConcentration(model, upstr);
        end
    end

end
 
