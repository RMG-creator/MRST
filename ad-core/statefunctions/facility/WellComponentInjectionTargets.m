classdef WellComponentInjectionTargets < StateFunction
    % Component total flux for wells (with treatment for cross-flow)
    properties

    end
    
    methods
        function gp = WellComponentInjectionTargets(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
%             gp = gp.dependsOn('ComponentPhaseFlux');
        end
        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            resmodel = model.ReservoirModel;
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            W   = map.W;
            ws = state.wellSol(map.active);
            targets = arrayfun(@(x) x.type, W, 'UniformOutput', false);
            val = vertcat(ws.val);
            nw = numel(W);
            surfaceComposition = cell(ncomp, nph);
            for c = 1:ncomp
                % Store well injector composition
                comp = resmodel.Components{c};
                surfaceComposition(c, :) = comp.getPhaseComponentFractionWell(resmodel, state, W);
            end
            rem = cellfun(@isempty, surfaceComposition);
            [surfaceComposition{rem}] = deal(zeros(nw, 1));
            compi = zeros(nw, ncomp);
            Wcomp = vertcat(W.compi);
            for ph = 1:nph
                compi = compi + Wcomp(:, ph).*[surfaceComposition{:, ph}];
            end
            
            isRateInjector = (strcmpi(targets, 'rate') & val >= 0);
            
            if isfield(W, 'rhoS')
                % Surface density is given on a per-well-basis for the
                % injectors
                rhoS = vertcat(W(isRateInjector).rhoS);
            else
                rhoS = resmodel.getSurfaceDensities();
                rhoS = repmat(rhoS(1, :), sum(isRateInjector), 1);
            end
            massRates = nan(nw, 1);
            massTargets = sum(bsxfun(@times, rhoS.*Wcomp(isRateInjector, :), val(isRateInjector)), 2);
            
            massRates(isRateInjector) = massTargets;
            v = struct('massRateTargets', massRates, 'surfaceMassComposition', compi, 'isRateInjector', isRateInjector);
        end
    end
end