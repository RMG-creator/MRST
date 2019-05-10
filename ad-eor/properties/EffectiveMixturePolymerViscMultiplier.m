classdef EffectiveMixturePolymerViscMultiplier < GridProperty
    properties
    end

    methods
        function gp = EffectiveMixturePolymerViscMultiplier(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'polymer'}, 'state'); % check mechanism
        end

        function muWeffMult = evaluateOnDomain(prop, model, state)
            c   = model.getProp(state, 'polymer');            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cbar   = c/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            b = 1./(1 - cbar + cbar./a);
            % The viscosity multiplier only results from the polymer mixing.
            muWeffMult = b.*fluid.muWMult(c).^mixpar;            
        end
    end
end