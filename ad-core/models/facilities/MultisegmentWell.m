classdef MultisegmentWell < SimpleWell
    properties
        signChangeChop
    end
    
    methods
        function well = MultisegmentWell(W, varargin)
            well = well@SimpleWell(W);
            assert(W.isMS)
            well.signChangeChop = false;
            if nargin > 1
                well = merge_options(well, varargin{:});
            end
            [nn, ns] = deal(numel(W.nodes.depth), numel(W.segments.length));
            C = sparse((1:ns)'*[1 1], W.segments.topo(1:end,:), ones(ns,1)*[1, -1], ns, nn);
            well.operators.grad = @(x)-C*x;
            well.operators.div  = @(x)C'*x;
            well.operators.segmentUpstr = @(flag, val)segmentUpstreamValue(flag, val, W.segments.topo);
            well.operators.C = C;
            averMat = sparse((1:ns)'*[1 1], W.segments.topo(1:end,:), ones(ns,1)*[1, ...
                                1], ns, nn);
            averMat = bsxfun(@rdivide, averMat, sum(averMat, 2));
            well.operators.aver = @(val) (averMat*val);
        end
        
        function counts = getVariableCounts(wm, fld)
            switch lower(fld)
                case {'pn', 'rw', 'ro', 'rg'}
                    counts = numel(wm.W.nodes.depth) - 1;
                case {'vmix'}
                    counts = size(wm.W.segments.topo, 1);
                otherwise
                    counts = getVariableCounts@SimpleWell(wm, fld);
            end
        end
        
        function [weqs, ctrlEq, weqsMS, extraNames, qMass, qSurf, wellSol] = computeWellEquations(well, wellSol0, wellSol, resmodel, q_s, bh, varw, pw, mobw, rhow, dissolved, compw, dt, iteration)
            % Node pressures for the well
            pN = varw{1};
            % Mass fraction for the phases
            alpha = varw(2:4);
            % Mixture mass flux in segments
            vmix = varw{5};
            
            % Create struct with the reservoir quantities
            resProps = struct();
            resProps.pressure = pw;
            resProps.mob = mobw;
            resProps.rho = rhow;
            resProps.dissolved = dissolved;
            resProps.components = compw;
            resProps.b = rhow;
            rhoS = resmodel.getSurfaceDensities();
            for i = 1:numel(resProps.b)
                % We store both b-factors and the density. This is a bit
                % redundant, but it makes it easier to write equations both
                % in terms of mass and surface volumes.
                den = rhoS(i);
                for j = 1:numel(resProps.b)
                    if ~isempty(dissolved{j}{i})
                        den = den + rhoS(j)*dissolved{j}{i};
                    end
                end
                resProps.b{i} = resProps.b{i}./den;
            end
            
            [weqs, weqsMS, qSurf, wellSol, mix_s, status, cstatus, qRes] = ...
                setupMSWellEquationSingleWell(well, resmodel, wellSol0, wellSol, q_s, bh, pN, alpha, vmix, resProps, dt, iteration);
            
            extraNames = well.getExtraEquationNames(resmodel);
            qMass = qSurf;
            for i = 1:numel(qMass)
                % Mass source terms
                qMass{i} = qSurf{i}.*rhoS(i);
            end
            ctrlEq =  setupWellControlEquationsSingleWell(wellSol, bh, q_s, status, mix_s, resmodel);
            % Update well properties which are not primary variables
            toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
            cq_sDb = cell2mat(toDouble(qSurf));
            
            wellSol.cqs     = cq_sDb;
            wellSol.cstatus = cstatus;
            wellSol.status  = status;
        end

        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = 1;
            switch(lower(name))
                case {'pn', 'nodepressure'}
                    fn = 'nodePressure';
                case {'rw', 'watermassfraction'}
                    fn = 'nodeComp';
                    index = 1;
                case {'ro', 'oilmassfraction'}
                    fn = 'nodeComp';
                    index = 2;
                case {'rg', 'gasmassfraction'}
                    fn = 'nodeComp';
                    index = 3;
                case {'vmix', 'segmentflux'}
                    fn = 'segmentFlux';
                    index = 1;
                otherwise
                    [fn, index] = getVariableField@SimpleWell(model, name);
            end
        end

        function ws = updateWellSol(well, ws, variables, dx)
            act = true(size(dx));
            for i = 1:numel(dx)
                name = variables{i};
                dv = dx{i};
                switch name
                    case 'pN'
                        dv = well.limitUpdateRelative(dv, ws.bhp, well.dpMaxRel);
                        dv = well.limitUpdateAbsolute(dv, well.dpMaxAbs);
                        ws.nodePressure = ws.nodePressure + dv;
                        fprintf('Biggest update %s: %f bar\n', name, max(abs(dv))/barsa);
                    case {'rW', 'rO', 'rG'}
                        ws = well.updateStateFromIncrement(ws, dv, [], name, inf, well.dsMaxAbs);
                        ws = well.capProperty(ws, name, 0, 1);
                        fprintf('Biggest update %s: %f\n', name, max(abs(dv)));
                    case 'vmix'
                        v = ws.segmentFlux;
                        if well.signChangeChop
                            v_changes_sign = ...
                                ((v + dv < 0) & (v >=0)) | ...
                                ((v + dv > 0) & (v <= 0));
                            if any(v_changes_sign)
                                v(v_changes_sign) = -eps*sign(v(v_changes_sign));
                                dv(v_changes_sign) = 0;
                            end
                        end
                        ws.segmentFlux = v + dv;
                        fprintf('Biggest relative update %s: %f\n', name, max(abs(dv))/norm(ws.segmentFlux, inf));
                    otherwise
                        continue
                end
                % We did update something
                act(i) = false;
            end
            ws = updateWellSol@SimpleWell(well, ws, variables(act), dx(act));
            ws = well.ensureWellSolConsistency(ws);
        end
        
        function names = getExtraPrimaryVariableNames(well, resmodel)
            names = getExtraPrimaryVariableNames@SimpleWell(well, resmodel);
            enames = {'rW', 'rO', 'rG'};
            enames = enames(resmodel.getActivePhases());
            names = [names, 'pN', enames, 'vmix'];
        end
        
        function [names, types] = getExtraEquationNames(well, resmodel)
            [names, types] = getExtraEquationNames@SimpleWell(well, resmodel);
            enames = {'waterNode', 'oilNode', 'gasNode'};
            enames = enames(resmodel.getActivePhases());
            names = [names, enames, 'pDropSeg', 'segMassClosure'];
            
            if nargout > 1
                etypes = {'node', 'node', 'node'};
                etypes = etypes(resmodel.getActivePhases());
                types = [types, etypes, 'seg', 'alpha'];
            end
        end
        
        function wellSol = validateWellSol(well, resmodel, wellSol, state)
            if ~isfield(wellSol, 'nodePressure') || isempty(wellSol.nodePressure)
                % Need to initialize the well
                
                W = well.W;
                nn  = numel(W.nodes.depth);
                %ns  = numel(W.segments.length);
                % set initial node pressures equal to initial bhp 
                wellSol.nodePressure = rldecode(wellSol.bhp, nn(:)-1);
                
                wc = W.cells;
                % get saturation and pressure from connecting cells
                sc = state.s(wc,:);
                pc = state.pressure(wc);
                % compute connecting cells props:
                [krw, kro, krg] = resmodel.relPermWOG(sc(:,1), sc(:,2), sc(:,3), resmodel.fluid);
                muw = resmodel.fluid.muW(pc);
                bw = resmodel.fluid.bW(pc);
                if resmodel.disgas
                    muo = resmodel.fluid.muO(pc, state.rs(wc), false(size(wc)));
                    bo  = resmodel.fluid.bO(pc, state.rs(wc), false(size(wc)));
                    rs  = state.rs(wc);
                else
                    muo = resmodel.fluid.muO(pc);
                    bo  = resmodel.fluid.bO(pc);
                    rs = 0;
                end
                if resmodel.vapoil
                    mug = resmodel.fluid.muG(pc, state.rv(wc), false(size(wc)));
                    bg  = resmodel.fluid.bG(pc, state.rv(wc), false(size(wc)));
                    rv  = state.rv(wc);
                else
                    mug = resmodel.fluid.muG(pc);
                    bg = resmodel.fluid.bG(pc);
                    rv = 0;
                end
                
                % phase mobilities
                [mw, mo, mg] = deal(krw./muw, kro./muo, krg./mug);

                % compute components weights for a unit volume rate  
                dens = resmodel.getSurfaceDensities();
                fwm = dens(1)*bw.*mw;
                fom = dens(2)*(bo.*mo + bg.*rv.*mg);
                fgm = dens(3)*(bg.*mg + bo.*rs.*mo);
                
                % devide by total to get component fractions in nodes
                % connected to reservoir grid cells
                nodemix = W.cell2node*bsxfun(@rdivide, [fwm, fom, fgm], fwm + fom + fgm) ;
       
                % set a small rate and compute a "plausible" mass conservative mixture rate
                % by solving least square problem 
                cq = -10*ones(numel(W.cells), 1)/day;
                cn = W.cell2node*cq;
                M = well.operators.C(:, 2:end);
                wellSol.segmentFlux =  -M*((M'*M)\cn(2:end));
                
                % set component mass fraction in undefined nodes equal to
                % average upstream node mass fractions
                nodemix_seg = -M*((M'*M)\nodemix(2:end,:));
                tp = W.segments.topo(:,2)-1;
                while any(sum(nodemix,2)==0)
                    nodemix(tp, :) = nodemix_seg;
                end
                nodemix = bsxfun(@rdivide, nodemix, sum(nodemix,2));
                nodemix = nodemix(2:end,:);
                
                wellSol.nodeComp = nodemix;
                % adjust total rates
                q_top = abs(well.operators.C(:,1)'*wellSol.segmentFlux);
                wellSol.qWs  = -q_top*nodemix(1,1)*dens(1);
                wellSol.qOs  = -q_top*nodemix(1,2)*dens(2);
                wellSol.qGs  = -q_top*nodemix(1,3)*dens(3);
                
                
            end
        end
        
        function ws = ensureWellSolConsistency(well, ws) %#ok
            % guarantees that the sum of rW, rO, rG remains equal to one.
            ws.nodeComp = bsxfun(@rdivide, ws.nodeComp, sum(ws.nodeComp, 2));
        end
        
    end
end

function v = segmentUpstreamValue(flag, val, topo)
    ix = flag.*topo(:,1) + ~flag.*topo(:,2);
    if ix(1) == 0
        ix(1) = topo(1,2);
    end
    v = val(max(ix-1,1));
end
