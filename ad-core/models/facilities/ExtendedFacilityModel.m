classdef ExtendedFacilityModel < FacilityModel
    properties
        T = 288.15; % Metric standard conditions
        pressure = 101.325*kilo*Pascal; % Metric standard pressure
        SeparatorGroup
    end
    
    methods
        function model = ExtendedFacilityModel(varargin)
            model@FacilityModel(varargin{:});
        end
        
        function n = getNumberOfComponents(fm)
            n = fm.ReservoirModel.getNumberOfComponents();
        end
        
        function n = getNumberOfPhases(fm)
            n = fm.ReservoirModel.getNumberOfPhases();
        end
        
        function src = getComponentSources(facility, state)
            map = facility.getProps(state, 'FacilityWellMapping');
            if isempty(map.W)
                val = [];
            else
                val = facility.getProps(state, 'ComponentTotalFlux');
            end
            src = struct('value', {val}, 'cells', map.cells);
        end
        
        function [surfaceRates, surfaceDensity] = getSurfaceRates(facility, state)
            [cflux, map] = facility.getProps(state, 'ComponentTotalFlux', 'FacilityWellMapping');
            model = facility.ReservoirModel;
            for c = 1:numel(cflux)
                % Sum over each well
                cflux{c} = map.perforationSum*cflux{c};
            end
            % We use a simple, but fast approach based on the
            % individual components' preference at different conditions
            [p, temp] = facility.getSurfaceConditions();
            nph = model.getNumberOfPhases();
            phaseMassRates = cell(1, nph);
            [phaseMassRates{:}] = deal(0);
            surfaceDensity = cell(1, nph);
            for c = 1:numel(cflux)
                composition = model.Components{c}.getPhaseCompositionSurface(model, state, p, temp);
                for ph = 1:nph
                    if ~isempty(composition{ph})
                        phaseMassRates{ph} = phaseMassRates{ph} + composition{ph}.*cflux{c};
                    end
                end
            end
            rhoS = facility.getProp(state, 'InjectionSurfaceDensity');
            surfaceRates = cell(1, nph);
            for ph = 1:nph
                rhoPhase = rhoS{ph};
                surfaceRates{ph} = phaseMassRates{ph}./rhoPhase;
                rhoPhase = model.AutoDiffBackend.convertToAD(rhoPhase, cflux{1});
                surfaceDensity{ph} = rhoPhase;
            end
            isProd = ~map.isInjector;
            if ~isempty(facility.SeparatorGroup) && any(isProd)
                % We have separators for the producers. Perform a
                % potentially costly separation to figure out phase rates.
                cfluxProd = cellfun(@(x) x(isProd), cflux, 'UniformOutput', false);
                [surfaceRatesProd, surfaceDensityProd] = facility.SeparatorGroup.getSurfaceRates(model, cfluxProd);
                for ph = 1:nph
                    surfaceRates{ph}(isProd) = surfaceRatesProd{ph};
                    surfaceDensity{ph}(isProd) = surfaceDensityProd{ph};
                end
            end
        end
        
        function [eqs, names, types, state] = getModelEquations(facility, state0, state, dt, drivingForces)
            model = facility.ReservoirModel;
            map = facility.getProps(state, 'FacilityWellMapping');
            if isempty(map.W)
                [eqs, names, types] = deal({});
                return
            end
            nph = model.getNumberOfPhases();
            wsum = map.perforationSum;
            % Get surface rate equations
            [surfaceRates, rhoS] = facility.getSurfaceRates(state);
            % Approximate scaling to get the tolerance in terms of
            % reservoir rates. Otherwise, e.g. gas phases may be very
            % differently scaled.
            rhoR = value(model.getProps(state0, 'Density'));
            rhoScale = bsxfun(@rdivide, value(rhoS), mean(rhoR, 1));
            % One equation for each phase corresponding to the volumetric
            % rate at surface conditions
            separateRates = strcmpi(facility.primaryVariableSet, 'standard');
            [sn, phnames] = model.getPhaseNames();
            if separateRates
                % This is a temporary hack!
                q_s = state.FacilityState.primaryVariables(1:nph);
                [eqs, names, types] = deal(cell(1, nph+1));
                for ph = 1:nph
                    eqs{ph} = (q_s{ph} - surfaceRates{ph}).*rhoScale(ph);
                    names{ph} = [phnames{ph}, 'Wells'];
                    types{ph} = 'perf';
                end
                targetRates = q_s;
                ctrl_index = nph+1;
            else
                varBHP = strcmpi(facility.primaryVariableSet, 'bhp');
                varMF = strcmpi(facility.primaryVariableSet, 'bhp_massfractions');
                assert(varBHP || varMF);
                % We need to actually store the surface rates in wellSol
                % here, since there are no corresponding primary variables
                if varMF
                    massFractions = state.FacilityState.massfractions;
                    cnames = model.getComponentNames();
                    componentSources = facility.getProps(state, 'ComponentTotalFlux');
                    ncomp = model.getNumberOfComponents();
                    [eqs, names, types] = deal(cell(1, ncomp));
                    total = 0;
                    for i = 1:ncomp
                        componentSources{i} = map.perforationSum*componentSources{i};
                        total = total + componentSources{i};
                    end
                    for i = 1:ncomp-1
                        eqs{i+1} = massFractions{i} - componentSources{i}./total;
                        names{i+1} = ['well_', cnames{i}];
                        types{i+1} = 'wellcomposition';
                    end
                else
                    [eqs, names, types] = deal(cell(1, 1));
                end
                ctrl_index = 1;
                targetRates = surfaceRates;
                qSurf = value(surfaceRates);
                for ph = 1:numel(sn)
                    fld = ['q', sn(ph), 's'];
                    for i = 1:numel(map.active)
                        ix = map.active(i);
                        state.wellSol(ix).(fld) = qSurf(i, ph);
                    end
                end
                clear qSurf;
            end
            % Set up AD for control equations
            nact = numel(map.active);
            bhp = state.FacilityState.primaryVariables{strcmpi(state.FacilityState.names, 'bhp')};
            backend = model.AutoDiffBackend;
            % Equation for well matching correct control
            ctrl_eq = backend.convertToAD(zeros(nact, 1), bhp);
            % Current control types (strings)
            well_controls = {state.wellSol(map.active).type}';
            % Current targets (numerical values)
            targets = vertcat(state.wellSol(map.active).val);

            % Handle bhp
            is_bhp = strcmp(well_controls, 'bhp');
            ctrl_eq(is_bhp) = bhp(is_bhp) - targets(is_bhp);

            % Following: Different types of rate controls.
            % Surface total rates
            is_rate = strcmp(well_controls, 'rate') | strcmpi(well_controls, 'vrat');
            % Surface oil rates
            is_orat = strcmp(well_controls, 'orat');
            % Surface water rates
            is_wrat = strcmp(well_controls, 'wrat');
            % Surface gas rates
            is_grat = strcmp(well_controls, 'grat');
            % Surface liquid rates (water + oil)
            is_lrat = strcmp(well_controls, 'lrat');
            % Reservoir rates (at averaged conditions from previous step)
            is_resv = strcmp(well_controls, 'resv');
            % Reservoir rates (at current conditions for each perf.)
            is_volume = strcmp(well_controls, 'volume');

            phases = model.getPhaseNames();
            is_surface_control = false(nact, 1);
            wrates = backend.convertToAD(zeros(nact, 1), bhp);
            for i = 1:nph
                switch phases(i)
                    case 'W'
                        act = is_rate | is_wrat | is_lrat;
                    case 'O'
                        act = is_rate | is_orat | is_lrat;
                    case 'G'
                        act = is_rate | is_grat;
                end
                is_surface_control(act) = true;
                wrates(act) = wrates(act) + targetRates{i}(act);
            end
            ctrl_eq(is_surface_control) = wrates(is_surface_control) - targets(is_surface_control);
            % RESV controls are special
            if any(is_resv)
                map = facility.getProp(state, 'FacilityWellMapping');
                rho = cellfun(@(x) x.ControlDensity, facility.WellModels(map.active), 'UniformOutput', false);
                rho = vertcat(rho{is_resv});
                resv_rates = 0;
                phaseRates = facility.getProps(state, 'PhaseFlux');
                rhoR = model.getProps(state, 'Density');
                for ph = 1:nph
                    tmp = wsum*(phaseRates{ph}.*rhoR{ph}(map.cells));
                    resv_rates = resv_rates + tmp(is_resv)./rho(:, ph);
                end
                ctrl_eq(is_resv) = resv_rates - targets(is_resv);
            end

            % Zero surface rate conditions
            wsign = vertcat(map.W.sign);
            % resv_value = value(resv_rates);
            surface_value = value(wrates);
            
            zeroTarget = targets == 0 & (is_surface_control | is_resv);
            % zeroSurface = (is_surface_control & sign(surface_value) ~= wsign & wsign ~= 0);
            % zeroRESV = (is_resv & sign(resv_value) ~= wsign & wsign ~= 0);
            zeroBHP = (is_bhp & sign(surface_value) ~= wsign & wsign ~= 0 & surface_value ~= 0);
%             zeroRates = zeroTarget | zeroSurface | zeroRESV | zeroBHP;
            zeroRates = zeroTarget | zeroBHP;
            if any(zeroRates)
                q_t = 0;
                for i = 1:nph
                    q_t = q_t + targetRates{i}(zeroRates);
                end
                ctrl_eq(zeroRates) = q_t;
            end

            if any(is_volume)
                phase_flux = facility.getProps(state, 'PhaseFlux');
                total_flux = 0;
                for i = 1:numel(phase_flux)
                    total_flux = total_flux + phase_flux{i};
                end
                well_total_flux = wsum*total_flux;
                ctrl_eq(is_volume) = well_total_flux(is_volume) - targets(is_volume);
            end

            assert(all(is_surface_control | is_bhp | is_volume | is_resv));

            eqs{ctrl_index} = ctrl_eq;
            names{ctrl_index} = 'closureWells';
            types{ctrl_index} = 'well';
        end
        
        function state = applyWellLimits(fm, state)
            active = fm.getIndicesOfActiveWells(state.wellSol);
            for i = 1:numel(active)
                w = active(i);
                well = fm.WellModels{w};
                state.wellSol(w) = fm.applyWellLimitsWellSol(well, state.wellSol(w));
            end
        end

        function model = validateModel(model, varargin)
            model = validateModel@FacilityModel(model, varargin{:});
        end

        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces)
            [model, state] = prepareReportstep@FacilityModel(model, state, state0, dt, drivingForces);
            [model, state] = model.updateRESVControls(state, state0);
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % Update pressure drop
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);
            if nw > 0
                rho = model.ReservoirModel.getProps(state, 'Density');
                rho = [rho{:}];
                [wc, perf2well] = model.getActiveWellCells(wellSol);
                rho = rho(wc, :);
                for i = 1:nw
                    wellNo = actWellIx(i);
                    wm = model.WellModels{wellNo};
                    % Possible index error for i here - should it be
                    % wellno?
                    rho_i = rho(perf2well == wellNo, :);
                    wellSol(wellNo) = wm.updateConnectionPressureDropState(model.ReservoirModel, wellSol(wellNo), rho_i, rho_i);
                end
            end
            state.wellSol = wellSol;
            state = model.applyWellLimits(state);
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Generic update function for reservoir models containing wells.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`

            [state, report] = updateAfterConvergence@FacilityModel(model, state0, state, dt, drivingForces);
            map = state.FacilityFluxProps.FacilityWellMapping;
            cf = state.FacilityFluxProps.ComponentTotalFlux;
            if iscell(cf)
                cf = [cf{:}];
            end
            for i = 1:numel(map.active)
                wi = map.active(i);
                act = map.perf2well == i;
                state.wellSol(wi).flux = state.FacilityFluxProps.PhaseFlux(act, :);
                state.wellSol(wi).ComponentTotalFlux = cf(act, :);
            end
        end
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            assert(not(isempty(model.FacilityFluxDiscretization)), ...
                'FacilityFluxDiscretization not initialized - did you call "validateModel"?');
            containers = [containers, {model.FacilityFluxDiscretization}];
        end
        
        function [p, T] = getSurfaceConditions(fm)
            p = fm.pressure;
            T = fm.T;
        end
        
        function wellSol = applyWellLimitsWellSol(fm, well, wellSol)
            % Update solution variables and wellSol based on the well
            % limits. If limits have been reached, this function will
            % attempt to re-initialize the values and change the controls
            % so that the next step keeps within the prescribed ranges.
            withinLimits = true;
            if ~well.allowControlSwitching
                % We cannot change controls, so we return
                return
            end
            if isfield(well.W, 'status') && ~well.W.status
                % Well is inactive
                return
            end
            lims = well.W.lims;
            model = fm.ReservoirModel;

            if ~isnumeric(lims)
                phases = model.getPhaseNames();
            
                qs_t = 0;
                for i = 1:numel(phases)
                    qs_t = qs_t + wellSol.(['q', phases(i), 's']);
                end
                bhp = wellSol.bhp;
            
                if well.isInjector()
                    % Injectors have three possible limits:
                    % bhp:  Upper limit on pressure.
                    % rate: Upper limit on total surface rate.
                    % vrat: Lower limit on total surface rate.
                    modes   = {'bhp', 'rate', 'vrat'};
                    lims = well.setMissingLimits(lims, modes, inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is lower limit, switch default sign
                        lims.vrat = -inf;
                    end

                    flags = [value(bhp)  > lims.bhp, ...
                              qs_t       > lims.rate, ...
                              qs_t       < lims.vrat];
                else
                    % Producers have several possible limits:
                    % bhp:  Lower limit on pressure.
                    % orat: Lower limit on surface oil rate
                    % lrat: Lower limit on surface liquid (water + oil) rate
                    % grat: Lower limit on surface gas rate
                    % wrat: Lower limit on surface water rate
                    % vrat: Upper limit on total volumetric surface rate

                    modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                    lims = fm.setMissingLimits(lims, modes, -inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is upper limit, switch default sign
                        lims.vrat = inf;
                    end
                    [q_w, q_o, q_g, q_sl] = deal(0);
                    if model.water
                        q_w = wellSol.qWs;
                    end
                    if model.oil
                        q_o = wellSol.qOs;
                    end
                    if model.gas
                        q_g = wellSol.qGs;
                    end
                    if isprop(model, 'solvent') && model.solvent
                        q_sl = wellSol.qSs;
                    end
                    flags = [value(bhp) < lims.bhp,  ...
                        q_o          < lims.orat, ...
                        q_w + q_o    < lims.lrat, ...
                        q_g + q_sl   < lims.grat, ...
                        q_w          < lims.wrat, ...
                        qs_t         > lims.vrat];
                end
            else
                modes = {};
                flags = false;
                assert(isempty(lims) || isinf(lims))
            end
            % limits we need to check (all others than w.type):
            chkInx = ~strcmp(wellSol.type, modes);
            vltInx = find(flags(chkInx), 1);
            if ~isempty(vltInx)
                withinLimits = false;
                modes  = modes(chkInx);
                switchMode = modes{vltInx};
                fprintf('Well %s: Control mode changed from %s to %s.\n', wellSol.name, wellSol.type, switchMode);
                wellSol.type = switchMode;
                wellSol.val  = lims.(switchMode);
            end

            if ~withinLimits
                v  = wellSol.val;
                switch wellSol.type
                    case 'bhp'
                        wellSol.bhp = bhp;
                    case 'rate'
                        for ix = 1:numel(phases)
                            wellSol.(['q', phases(ix), 's']) = v*well.W.compi(ix);
                        end
                    case 'orat'
                        wellSol.qOs = v;
                    case 'wrat'
                        wellSol.qWs = v;
                    case 'grat'
                        wellSol.qGs = v;
                end % No good guess for qOs, etc...
            end
        end
        
        function lims = setMissingLimits(fm, lims, modes, val)
            missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
            for f = missing_fields
               lims = setfield(lims, f{:}, val);
            end
        end
        
        function [model, state] = updateRESVControls(model, state, state0)
            % Treat RESV
            activeWellMask = model.getWellStatusMask(state.wellSol);
            isRESVHist = cellfun(@(x) strcmpi(x.W.type, 'resv_history'), model.WellModels(activeWellMask));
            isRESVNow = cellfun(@(x) strcmpi(x.W.type, 'resv'), model.WellModels(activeWellMask));
            isRESV = isRESVHist | isRESVNow;
            if any(isRESV)
                % Local index
                W = model.getWellStruct(activeWellMask);
                
                isHist = isRESVHist(isRESV);
                compi = vertcat(W.compi);
                compi = compi(isRESV, :);
                
                rates = vertcat(W(isRESV).val);
                qs = bsxfun(@times, rates, compi);
                
                rmodel = model.ReservoirModel;
                disgas = isprop(rmodel, 'disgas') && rmodel.disgas;
                vapoil = isprop(rmodel, 'vapoil') && rmodel.vapoil;
                oix = rmodel.getPhaseIndex('O');
                gix = rmodel.getPhaseIndex('G');
                
                pvt_reg = rmodel.FlowPropertyFunctions.Density.regions;
                if isempty(pvt_reg)
                    pvt_reg = ones(rmodel.G.cells.num, 1);
                end
                cells = arrayfun(@(x) x.cells(1), W(isRESV));
                nc = numel(cells);
                regNo = pvt_reg(cells);
                
                rs = zeros(nc, 1);
                rv = zeros(nc, 1);
                pm = zeros(nc, 1);
                pv = rmodel.getProp(state0, 'PoreVolume');
                if rmodel.water
                    sw = rmodel.getProp(state0, 'sw');
                    pv = pv.*(1-sw);
                end
                for reg = 1:max(regNo)
                    subs = regNo == reg;
                    local = pvt_reg == reg;
                    pvi = pv.*local;
                    pm(subs) = sum(state0.pressure.*pvi)/sum(pvi);
                    if disgas
                        pvi = pv.*local;
                        rs(subs) = sum(state0.rs.*pvi)/sum(pvi);
                    end
                    if vapoil
                        pvi = pv.*local;
                        rv(subs) = sum(state0.rv.*pvi)/sum(pvi);
                    end
                end
                if any(isHist)
                    if disgas
                        rs(isHist) = min(qs(isHist, gix)./qs(isHist, oix), rs);
                    end
                    if vapoil
                        rv(isHist) = min(qs(isHist, oix)./qs(isHist, gix), rv);
                    end
                end
                substate = struct('pressure', pm, ...
                                  's', repmat([1, 0, 0], nc, 1), ...
                                  'rs', rs, ...
                                  'rv', rv);

                flowProps = model.ReservoirModel.FlowPropertyFunctions.subset(cells);
                % Avoid using flag for interpolation
                flowProps.ShrinkageFactors.useSaturatedFlag = true;
                substate = flowProps.evaluateStateFunctionWithDependencies(model.ReservoirModel, substate, 'Density');
                rho = substate.FlowProps.Density;
                rho = [rho{:}];
                if false
                    newRates = sum(qs.*rhoS./rho, 2);
                else
                    shrink = 1 - rs.*rv;
                    b = substate.FlowProps.ShrinkageFactors;
                    newRates = 0;
                    if rmodel.water
                        wix = 1;
                        bW = b{wix};
                        newRates = newRates + qs(:, wix)./bW;
                    end
                    if rmodel.oil
                        bO = b{oix};
                        orat = qs(:, oix);
                        if vapoil
                            orat = orat - rv.*qs(:, gix);
                        end
                        newRates = newRates + orat./(bO.*shrink);
                    end
                    if rmodel.gas
                        bG = b{gix};
                        grat = qs(:, gix);
                        if vapoil
                            grat = grat - rs.*qs(:, oix);
                        end
                        newRates = newRates + grat./(bG.*shrink);
                    end
                end
                resvIx = find(isRESV);
                actIx = find(activeWellMask);
                for i = 1:numel(resvIx)
                    I = resvIx(i);
                    global_well_ix = actIx(I);
                    if isRESVHist(I)
                        model.WellModels{global_well_ix}.W.val = newRates(i);
                        state.wellSol(global_well_ix).val = newRates(i);
                        model.WellModels{global_well_ix}.W.type = 'resv';
                        state.wellSol(global_well_ix).type = 'resv';
                    end
                    model.WellModels{global_well_ix}.ControlDensity = rho(i, :);
                end
            end
        end
    end
end