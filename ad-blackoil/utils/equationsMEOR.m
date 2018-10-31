function [problem, state] = equationsMEOR(state0,state,model,dt,drivingForces,varargin)
% Work in progress to create MEOR effects
% Get linearized problem for oil/water/MEOR system with black oil
% properties
opt = struct('Verbose', mrstVerbose,...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);
         
opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
%assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

% Operators, grid, and fluid model
s = model.operators;
G = model.G;
f = model.fluid;

% Properties at current timestep
[p, sW, m, n, wellSol] = model.getProps(state, 'pressure', 'water',...
    'microbe', 'nutrient', 'wellsol');

% Properties at previous timestep
[p0,sW0, m0, n0] = model.getProps(state0, 'pressure', 'water', ...
    'microbe', 'nutrient');

pBH = vertcat(wellSol.bhp);
qWs = vertcat(wellSol.qWs);
qOs = vertcat(wellSol.qOs);
qWMEOR = vertcat(wellSol.qWMEOR);
%qWMICROBE = vertcat(wellSol.qWMICROBE);
%qWNUTRIENT = vertcat(wellSol.qWNUTRIENT);

% Initialize independent variables
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, m, n, qWs, qOs, qWMEOR, pBH] = ...
            initVariablesADI( p, sW, m, n, qWs, qOs, qWMEOR, pBH);
    else
        [p0, sW0, m0, n0, tmp,tmp,tmp,tmp] = ...
            initVariablesADI(p0, sW0, m0, n0,...
            zeros(size(qWs)), zeros(size(qOs)), zeros(size(qWMEOR)),...
            zeros(size(pBH)));
        clear tmp
    end
end

% We will solve for pressure, water saturation (oil saturation follows from
% the definition of saturations) ((may need to change later)), microbe
% concentration, nutrient concentration, later metabolite concentration,
% and well rates and bhp.
primaryVars = {'pressure', 'sW', 'microbe', 'nutrient', 'qWs', 'qOs',...
    'qWMEOR', 'bhp'};

%% old code
% % Pressure dependent transmissibility multiplier
% [trMult, pvMult, pvMult0, transMult] = deal(1);
% if isfield(f, 'tranMultR')
%     trMult = f.tranMultR(p);
% end
% % Pressure dependent pore volume multiplier (add biofilm here?)
% if isfield(f,'pvMultR')
%     pvMult = f.pvMultR(p);
%     pvMult0 = f.pvMultR(p0);
% end
% if isfield(f, 'transMult')
%     transMult = f.transMult(p);
% end
% % Check for capillary pressure (p_cow) (do relperm here?)
% pcOW = 0;
% if isfield(f, 'pcOW')
%     pcOW = f.pcOW(sW);
% end
% % Gravity contribution, assert that it is aligned with z-direction
% grav = gravity();
% assert(grav(1) == 0 && grav(2) == 0);
% g = norm(grav);
% dz = s.grad(G.cells.centroids(:,3));

%% More old code
%bW = f.bW(p);
%rhoW = bW.*f.rhoWS;
% Defining rhoW on faces as average of neighboring cells
%rhoWf = s.faceAvg(rhoW);
%muW = f.muW(p-pcOW);
%mobW = trMult.*krW./muW;
%dpW = s.grad(p-pcOw) - g*(rhoWf.*dz);
% Water upstream index
%upcw = (double(dpW)>=0);
%vW = s.faceUpstr(upcw, mobW).*T.*dpW;

%if any(bW<0)
    %warning('Negative water compressibility present')
%end

% Microbes and Nutrients (will change m to metabolites)
%mobM = (mobW.*m);
%mobN = (mobW.*n);
%vM = s.faceUpstr(upcw, mobM).*s.T.*dpW;
%vN = s.faceUpstr(upcw, mobN).*s.T.*dpW;

% Oil props
%bO = f.bO(p);
%rhoO = bO.*f.rhoOS;
%rhoOf = s.faceAvg(rhoO);
%dpO = s.grad(p) - g*(rhoOf.*dz);
% Oil upstream index
%upco = (double(dpO)>=0);
%if isfield(f, 'BOxmuO')
    % mob0 is already multipied with b0
 %   mobO = trMult.*krO./f.BOxmuO(p);
  %  bOvO = s.faceUpstr(upco, mobO).*T.*dpO;
  %  vO = bOvO./s.faceUpstr(upco, bO);
%else
 %   mobO = trMult.*krO./f.muO(p);
 %   vO = s.faceUpstr(upco, mobO).*T.*dpO;

%end
%if any(bO<0)
 %   warning('Negative oil compressibility present')
%end

%%
% Evaluate relative permeability
sO = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO}); % as written

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p,p0);

% Modify relperm by mobility multiplier
krW = mobMult.*krW;
krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water and MEOR props
[vW, vM, vN, bW, mobW, mobM,mobN, rhoW, pW, upcw] = ...
    getFluxAndPropsMEOR(model, p, sW, m, n, krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate Oil properties
[vO,bO,mobO,rhoO,p,upco] = getFluxAndPropsOil_BO(model,p,sO,krO,T,gdz);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vM, vN); %add function to model
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobM, mobN);%add function
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end


% EQUATIONS ---------------------------------------------------------------

bWvW = s.faceUpstr(upcw, bW).*vW;
bWvM = s.faceUpstr(upcw, bW).*vM;
bWvN = s.faceUpstr(upcw, bW).*vN;
bOvO = s.faceUpstr(upco, bO).*vO;


% Conservation of oil:
oil = (s.pv/dt).*(pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);
        
%Conservation of  water:
water = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*bW0.*sW0) + s.Div(bWvW);
        
%Conservation of  microbes:
microbe = (s.pv/dt).*(pvMult.*bW.*sW.*m - pvMult0.*bW0.*sW0.*m0)+ s.Div(bWvM);
        
%Conservation of  nutrients:
nutrient = (s.pv/dt).*(pvMult.*bW.*sW.*n - pvMult0.*bW0.*sW0.*n0)+ s.Div(bWvN);


eqs = {water, oil, microbe, nutrient}; 
names = {'water', 'oil', 'microbe', 'nutrient'};
types = {'cell', 'cell', 'cell', 'cell'};

% Add in any fluxes/source terms given as boundary conditions
[eqs, qBC, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
    model, eqs, {pW, p}, {rhoW, rhoO}, {mobW, mobO}, {bW, bO}, ...
    {sW, sO}, drivingForces);

% Add MEOR boundary conditions
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc,'microbe')
    injInx = qBC{1} > 0; % Water inflow indicies
    mbc = (BCTocellMap')*m;
    nbc = (BCTocellMap')*n;
    mbc(injInx) = drivingForces.bc.m(injInx);
    nbc(injInx) = drivingForces.bc.n(injInx);
    eqs{3} = eqs{3} - BCTocellMap*(mbc.*qBC{1});
    eqs{4} = eqs{4} - BCTocellMap*(nbc.*qBC{1});
end

% Add MEOR source
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'm')
    injInx = qSRC{1}>0;
    msrc = m(srcCells);
    nsrc = n(srcCells);
    msrc(injInx) = drivingForces.src.m(injInx);
    nsrc(injInx) = drivingForces.src.n(injInx);
    eqs{3}(srcCells) = eqs{3}(srcCells) - msrc.*qSRC{1};
    eqs{4}(srcCells) = eqs{4}(srcCells) - nsrc.*qSRC{1};
end

% well equations switch to 8 eqns with only 1 meor well
if ~isempty(W)
    wm = model.wellmodel;
    if ~opt.reverseMode
        wc = vertcat(W.cells);
        pw = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw = {bW(wc), bO(wc)};
        tw = {mobW(wc), mobO(wc)};
        s = {sW(wc), sO(wc)};
        
        [cqs, weqs, ctrleqs, wc, state.wellSol] = ...
            wm.computeWellFlux(model, W, wellSol,...
            pBH, {qWs, qOs}, pw, rhos, bw, tw, s, {},...
            'nonlinearIteration', opt.iteration);
        
        % Store the well equations (relating well BHP to influx)
        eqs(5:6) = weqs;
        
        % Store control equations
        eqs{8} = ctrleqs;
        % Add source terms to the equations.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{1}(wc) - cqs{2};
        
        % MEOR well equations
         [~, wciMEOR, iInxW, MEORc] = getWellMEOR(W);
         mw = m(wc);
         nw = n(wc);
         mw(iInxW) = wciMEOR.*(1-MEORc);
         nw(iInxW) = wciMEOR.*MEORc;
%        [~, wciMICROBE, iInxWM] = getWellMicrobe(W);
%        [~, wciNUTRIENT, iInxWN] = getWellNutrient(W);
%        mw = m(wc);
%        nw = n(wc);
%        mw(iInxWM) = wciMICROBE;
%        nw(iInxWN) = wciNUTRIENT;
        
        % Totally unsure here
        bWqM = mw.*cqs{1};
        bWqN = nw.*cqs{1};
        eqs{3}(wc) = eqs{3}(wc) - bWqM;
        eqs{4}(wc) = eqs{4}(wc) - bWqN;
        
        % Well MEOR rate for each well is water rate in each perforation
        % multiplied with microbe and nutrient concentration in that
        % perforated cell
        perf2well = getPerforationToWellMapping(W);
        Rw = sparse(perf2well, (1:numel(perf2well))', 1,numel(W),numel(perf2well));
        eqs{7} = qWMEOR - Rw*(cqs{1}.*(mw+nw)); %UNSURE
%        eqs{7} = qWMICROBE - Rw*(cqs{1}.*mw);
%        eqs{8} = qWNUTRIENT - Rw*(cqs{1}.*nw);
        
        names(5:8) = {'waterWells', 'oilWells', 'meorWells', 'closureWells'};
        types(5:8) = {'perf', 'perf', 'perf', 'well'};
    else
        [eq, n, typ] = ...
            wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        % add another equation for MEOR well rates. No idea if functional
        [eqs{5:8}] = deal(eq{1});
        [names{5:8}] = deal(n{1});
        [types{5:8}] = deal(typ{1});
    end
end
problem = LinearizedProblem(eqs,types,names,primaryVars,state,dt);
end

function [wMICROBE, wciMICROBE, iInxWM] = getWellMicrobe(W)
    if isempty(W)
        wMICROBE = [];
        wciMICROBE = [];
        iInxWM = [];
        return
    end
    inj = vertcat(W.sign) == 1;
    mInj = cellfun(@(x)~isempty(x), {W(inj).MICROBE});
    wMICROBE = zeros(nnz(inj),1);
    wMICROBE(mInj) = vertcat(W(inj(mInj)).MICROBE);
    wciMICROBE = rldecode(wMICROBE, cellfun(@numel, {W(inj).cells}));
    
    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw = numel(W);
    perf2well = rldecode((1:nw)',nPerf);
    compi = vertcat(W.compi);
    iInx = rldecode(inj, nPerf);
    iInx = find(iInx);
    iInxWM = iInx(compi(perf2well(iInx),1)==1);
end

function [wNUTRIENT, wciNUTRIENT, iInxWN] = getWellNutrient(W)
    if isempty(W)
        wNUTRIENT = [];
        wciNUTRIENT = [];
        iInxWN = [];
        return
    end
    inj = vertcat(W.sign) == 1;
    mInj = cellfun(@(x)~isempty(x), {W(inj).MICROBE});
    wNUTRIENT = zeros(nnz(inj),1);
    wNUTRIENT(mInj) = vertcat(W(inj(mInj)).MICROBE);
    wciNUTRIENT = rldecode(wNUTRIENT, cellfun(@numel, {W(inj).cells}));
    
    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw = numel(W);
    perf2well = rldecode((1:nw)',nPerf);
    compi = vertcat(W.compi);
    iInx = rldecode(inj, nPerf);
    iInx = find(iInx);
    iInxWN = iInx(compi(perf2well(iInx),1)==1);
end

function [wMEOR, wciMEOR, iInxW,MEORc] = getWellMEOR(W)
    if isempty(W)
        wMEOR = [];
        wciMEOR = [];
        iInxW = [];
        MEORc = [];
        return
    end
    inj = vertcat(W.sign) == 1;
    mInj = cellfun(@(x)~isempty(x), {W(inj).meor});
    wMEOR = zeros(nnz(inj),1);
    wMEOR(mInj) = vertcat(W(inj(mInj)).meor);
    wciMEOR = rldecode(wMEOR, cellfun(@numel, {W(inj).cells}));
    MEORcomp = zeros(nnz(inj),1);
    %MEORcomp(mInj) = vertcat(W(inj(mInj)).MEORMIX);
    MEORcomp(mInj) = [.5];
    MEORc = rldecode(MEORcomp, cellfun(@numel, {W(inj).cells}));
    
    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw = numel(W);
    perf2well = rldecode((1:nw)',nPerf);
    compi = vertcat(W.compi);
    iInx = rldecode(inj, nPerf);
    iInx = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end
        
        