classdef PolymerViscMultiplier < GridProperty
    properties
    end

    methods
        function gp = PolymerViscMultiplier(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'PolymerAdsorption'});
            gp = gp.dependsOn({'polymer'}, 'state'); % check mechanism
        end

        function muWMult = evaluateOnDomain(prop, model, state)

            c   = model.getProp(state, 'polymer');
            ads = model.getProp(state, 'PolymerAdsorption');
            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cbar   = c/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            b = 1./(1 - cbar + cbar./a);
            % The viscosity multiplier only results from the polymer mixing.
            muWeffMult = b.*fluid.muWMult(c).^mixpar;

            permRed = 1 + ((fluid.rrf - 1)./fluid.adsMax).*ads;
            muWMult  = muWeffMult.*permRed;
            
        end
    end
end