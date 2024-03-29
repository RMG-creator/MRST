classdef SurfactantAdsorption < StateFunction
    properties
    end

    methods
        function prop = SurfactantAdsorption(model, varargin)
            prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'surfactant', 'surfactantmax'}, 'state'); % check mechanism
        end

        function ads = evaluateOnDomain(prop, model, state)
            [cs, csmax] = model.getProps(state, 'surfactant', 'surfactantmax');
            fluid = model.fluid;
            ads  = effads(cs, csmax, fluid);
        end
    end
end

function y = effads(cs, csmax, f)
   if f.adsInxSft == 2
      y = f.surfads(max(cs, csmax));
   else
      y = f.surfads(cs);
   end
end