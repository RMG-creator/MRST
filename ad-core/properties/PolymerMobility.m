classdef PolymerMobility < GridProperty
    properties
    end
    
    methods
        function gp = Mobility(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn({'RelativePermeability', 'Viscosity'});
            gp = gp.dependsOn({'polymer'}, 'state'); % check mechanism
        end
        function mob = evaluateOnDomain(prop, model, state)
            [mu, kr] = prop.getEvaluatedDependencies(state, 'Viscosity', 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            
            
            %% polymer part
            fluid = model.fluid;
            c = model.getProp(state, 'polymer');
            mixpar = fluid.mixPar;
            cbar   = c/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            b = 1./(1-cbar+cbar./a);
            % The viscosity multiplier only results from the polymer mixing.
            muWeffMult = b.*fluid.muWMult(c).^mixpar;

            % TODO later: add adsorption effect
            % permRed = 1 + ((fluid.rrf-1)./fluid.adsMax).*ads;
            % muWMult  = muWeffMult.*permRed;
            
            muWMult  = muWeffMult.*permRed;
            % assuming water index is 1
            mob{1} = mob{1}/muWMult; % get water index in the proper way
            %% End polymer part 
            
            
            
            if isfield(model.fluid, 'tranMultR')
                % Pressure dependent mobility multiplier 
                p = model.getProp(state, 'pressure');
                mult = model.fluid.tranMultR(p);
                mob = cellfun(@(x) x.*mult, mob, 'UniformOutput', false);
            end
            mv = cellfun(@(x) min(value(x)), mob);
            if any(mv < 0)
            	warning('Negative mobilities detected! Capping to zero.')
                mob = cellfun(@(x) max(x, 0), mob, 'UniformOutput', false);
            end
        end
    end
end