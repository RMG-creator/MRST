classdef PolymerPhaseFlux < GridProperty
    properties

    end
    
    methods
        function gp = PolymerPhaseFlux(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp = gp.dependsOn({'FaceMobility', 'FaceConcentration', 'PolymerViscMultiplier', 'PhaseFlux'});
        end
        
        function v = evaluateOnDomain(prop, model, state)
          
            [mob, cf, muWMult, phaseFlux] = prop.getEvaluatedDependencies(state,...
                'FaceMobility', 'FaceConcentration','PolymerViscMultiplier', 'PhaseFlux');
            nph = numel(mob);
            v = cell(1, nph+1);
            for i = 1:nph
                v{i} = deal(phaseFlux{i});
            end
            
            v{1} = v{1}./muWMult;
            c   = model.getProp(state, 'polymer');
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cbar   = c/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            v{nph+1} = v{1}./(1+(1-a)*cbar).*cf;
        end
    end
end