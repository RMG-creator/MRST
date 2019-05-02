classdef PolymerAdsorption < GridProperty
    properties
    end

    methods
        function gp = PolymerAdsorption(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'polymer', 'polymermax'}, 'state'); % check mechanism
        end

        function ads = evaluateOnDomain(prop, model, state)
            [c, cmax] = model.getProps(state, 'polymer', 'polymermax');
            fluid = model.fluid;
            ads  = effads(c, cmax, fluid);
        end
    end
end

function y = effads(c, cmax, f)
   if f.adsInx == 2
      y = f.ads(max(c, cmax));
   else
      y = f.ads(c);
   end
end