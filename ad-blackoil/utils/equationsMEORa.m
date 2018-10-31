function [problem, state] = equationsMEORa(state0,state,model,dt,drivingForces,varargin)
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
Y_micro = model.yield_microbe;
Y_meta = model.yield_metabolite;
mu_micro = model.growth_max_microbe;
mu_meta = model.growth_max_metabolite;
K_micro = model.halfsat_microbe;
K_meta = model.halfsat_metabolite;
N = model.crit_val;


% Properties at current timestep
[p, sW, m, meta, n, wellSol] = model.getProps(state, 'pressure', 'water',...
    'microbe', 'metabolite', 'nutrient', 'wellsol');

% Properties at previous timestep
[p0,sW0, m0, meta0, n0] = model.getProps(state0, 'pressure', 'water', ...
    'microbe', 'metabolite', 'nutrient');

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
        [p, sW, m, meta, n, qWs, qOs, qWMEOR, pBH] = ...
            initVariablesADI( p, sW, m, meta, n, qWs, qOs, qWMEOR, pBH);
    else
        [p0, sW0, m0, meta0, n0, tmp,tmp,tmp,tmp] = ...
            initVariablesADI(p0, sW0, m0, meta0, n0,...
            zeros(size(qWs)), zeros(size(qOs)), zeros(size(qWMEOR)),...
            zeros(size(pBH)));
        clear tmp
    end
end

% We will solve for pressure, water saturation (oil saturation follows from
% the definition of saturations) ((may need to change later)), microbe
% concentration, nutrient concentration, metabolite concentration,
% and well rates and bhp.
primaryVars = {'pressure', 'sW', 'microbe', 'metabolite', 'nutrient', 'qWs', 'qOs',...
    'qWMEOR', 'bhp'};

% Evaluate relative permeability
sO = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO}); % as written

if model.biosurf
    % form is taken from Nielsen
    % constants are as well
    partition = 1.*(sW.*model.fluid.rhoWS)./(sO.*model.fluid.rhoOS);
    meta_eff = meta.*partition./(partition + 1);
    %meta_eff = meta;
    surfa = [1*10^-4, 0.2, 1.5*10^4];
    ift = @(s) 29.*(-tanh(surfa(3).*s-surfa(2))+1+surfa(1))./...
        (-tanh(-surfa(2))+1+surfa(1));
    sigma = ift(double(meta_eff));
    f = (sigma/29).^(1/6);
    sor = f.*.4;
    %disp(min(sor))
    f = ones(length(f),1);
    wmax = f.*.5+1-f;
    swi = f.*.3;
    omax = f.*.8+1-f;
    a = f.*2+1-f;
    krO_max = omax.*((sO-sor)./(1-swi-sor)).^a;
    krW_max = wmax.*((sW-swi)./(1-swi-sor)).^a;
    inx = meta_eff>1e-16;
    
%     krO_max =  .8.*((1 - sW - .08)./(1 - .08 - .3)).^2;
%     krW_max =  .5.*((sW - .3)./(1 - .08 - .3)).^2;
%     meta_eff = meta.*10^3;
%     partition = 1.*(sW.*model.fluid.rhoWS)./(sO.*model.fluid.rhoOS);
%     meta_eff = meta_eff.*partition./(partition + 1);
%     meta_eff(meta_eff<.01) = 0;
%     meta_eff(meta_eff>1) = 1;
     krW = krW + (krW_max - krW).*inx;
     krO = krO + (krO_max - krO).*inx;
    %krW = krW + (sW - krW).*meta_eff;
    %krO = krO + (sO - krO).*meta_eff;
end

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
[vW, vMicro, vMeta, vN, bW, mobW, mobM_0, mobPol, mobN, rhoW, pW, upcw] = ...
    getFluxAndPropsMEORa(model, p, sW, m, meta, n, krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate Oil properties
[vO,bO,mobO,rhoO,p,upco] = getFluxAndPropsOil_BO(model,p,sO,krO,T,gdz);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vMicro); %add function to model for more than 1 extra phase
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobM_0, mobN);%add function if I feel like it same as storeFluxes
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end


% EQUATIONS ---------------------------------------------------------------
% Microbe and metabolite calculations (should they be explicit or in eqn?)
mu_b = mu_micro.*n./(K_micro + n);
mu_m = mu_meta.*(n-N)./(K_meta + n - N);
R_n = -mu_b.*m.*Y_micro - mu_m.*m.*Y_meta;



bWvW = s.faceUpstr(upcw, bW).*vW;
bWvMicro = s.faceUpstr(upcw, bW).*vMicro;
bWvMeta = s.faceUpstr(upcw, bW).*vMeta;
bWvN = s.faceUpstr(upcw, bW).*vN;
bOvO = s.faceUpstr(upco, bO).*vO;


% Conservation of oil:
oil = (s.pv/dt).*(pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);
        
%Conservation of  water:
water = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*bW0.*sW0) + s.Div(bWvW);
        
%Conservation of  microbes:
microbe = (s.pv/dt).*((pvMult.*bW.*sW.*m - pvMult0.*bW0.*sW0.*m0)) - s.pv.*mu_b.*m.*bW.*sW.*pvMult.*Y_micro + s.Div(bWvMicro);

%Conservation of  nutrients:
nutrient = (s.pv/dt).*((pvMult.*bW.*sW.*n - pvMult0.*bW0.*sW0.*n0)) - s.pv.*R_n.*bW.*sW.*pvMult+ s.Div(bWvN);

%Conservation of metabolites:
metabolite = (s.pv/dt).*((pvMult.*bW.*sW.*meta - pvMult0.*bW0.*sW0.*meta0)) - s.pv.*mu_m.*m.*bW.*sW.*pvMult.*Y_meta+ s.Div(bWvMeta);

eqs = {water, oil, microbe, nutrient, metabolite}; 
names = {'water', 'oil', 'microbe', 'nutrient', 'metabolite'};
types = {'cell', 'cell', 'cell', 'cell', 'cell'};


% Add microbe processes before boundary conditions? or after?
% before equations? as equations?


% Add in any fluxes/source terms given as boundary conditions
[eqs, qBC, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(...
    model, eqs, {pW, p}, {rhoW, rhoO}, {mobW, mobO}, {bW, bO}, ...
    {sW, sO}, drivingForces);

% Add MEOR boundary conditions
if ~isempty(drivingForces.bc) && isfield(drivingForces.bc,'m')
    injInx = qBC{1} > 0; % Water inflow indicies
    mbc = (BCTocellMap')*m; % m_0 is only type injected
    nbc = (BCTocellMap')*n;
    metabc = (BCTocellMap')*meta;
    mbc(injInx) = drivingForces.bc.m(injInx);
    nbc(injInx) = drivingForces.bc.n(injInx);
    metabc(injInx) = drivingForces.bc.meta(injInx);
    eqs{3} = eqs{3} - BCTocellMap*(mbc.*qBC{1});
    eqs{4} = eqs{4} - BCTocellMap*(nbc.*qBC{1});
    eqs{5} = eqs{5} - BCTocellMap*(metabc.*qBC{1});
end

% Add MEOR source
if ~isempty(drivingForces.src) && isfield(drivingForces.src, 'm')
    injInx = qSRC{1}>0;
    msrc = m(srcCells);
    nsrc = n(srcCells);
    metasrc = meta(srcCells);
    msrc(injInx) = drivingForces.src.m(injInx);
    nsrc(injInx) = drivingForces.src.n(injInx);
    eqs{3}(srcCells) = eqs{3}(srcCells) - msrc.*qSRC{1};
    eqs{4}(srcCells) = eqs{4}(srcCells) - nsrc.*qSRC{1};
 %   eqs{5}(srcCells) = eqs{5}(srcCells) - metasrc.*qSRC{1};
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
        eqs(6:7) = weqs;
        
        % Store control equations
        eqs{9} = ctrleqs;
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
        eqs{4}(wc) = eqs{5}(wc) - bWqN;
        
        % Well MEOR rate for each well is water rate in each perforation
        % multiplied with microbe and nutrient concentration in that
        % perforated cell
        perf2well = getPerforationToWellMapping(W);
        Rw = sparse(perf2well, (1:numel(perf2well))', 1,numel(W),numel(perf2well));
        eqs{8} = qWMEOR - Rw*(cqs{1}.*(mw+nw)); %UNSURE
%        eqs{7} = qWMICROBE - Rw*(cqs{1}.*mw);
%        eqs{8} = qWNUTRIENT - Rw*(cqs{1}.*nw);
        
        names(6:9) = {'waterWells', 'oilWells', 'meorWells', 'closureWells'};
        types(6:9) = {'perf', 'perf', 'perf', 'well'};
    else
        [eq, n, typ] = ...
            wm.createReverseModeWellEquations(model, state0.wellSol, p0);
        % add another equation for MEOR well rates. No idea if functional
        [eqs{6:9}] = deal(eq{1});
        [names{6:9}] = deal(n{1});
        [types{6:9}] = deal(typ{1});
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
        
        