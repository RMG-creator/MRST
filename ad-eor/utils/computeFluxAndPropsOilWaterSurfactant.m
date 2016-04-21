function [dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T] = ...
        computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, pBH, W, varargin)
    
    opt = struct('velocCompMethod', 'square');
    opt = merge_options(opt, varargin{:});

    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;
    
    Nc = computeCapillaryNumber(p, c, pBH, W, fluid, G, op, 'velocCompMethod', opt.velocCompMethod);
    [krW, krO] = computeRelPermSft(sW, c, Nc, fluid);

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Compute transmissibility
    T = op.T.*transMult;

    % Gravity contribution
    gdz = model.getGravityGradient();

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW;
    krO = mobMult.*krO;

    fluid = model.fluid;

    % Capillary pressure
    pcOW = 0;
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
    end
    pcOW = pcOW.*fluid.ift(c)/fluid.ift(0);
    pO = p;
    pW = pO - pcOW;

    bW0 = fluid.bW(p0);
    bO0 = fluid.bO(p0);

    % Water flux and properties
    bW      = fluid.bW(p);
    rhoW    = bW.*fluid.rhoWS;
    rhoWf   = op.faceAvg(rhoW);
    muW     = fluid.muWSft(c);
    multmuW = fluid.muW(p)/fluid.muWr;
    mobW    = krW./(muW.*multmuW);
    dpW     = op.Grad(pW) - rhoWf.*gdz;
    upcw    = (double(dpW)<=0);

    % Oil flux and properties
    bO    = fluid.bO(pO);
    rhoO  = bO.*fluid.rhoOS;
    rhoOf = op.faceAvg(rhoO);
    if isfield(fluid, 'BOxmuO')
        muO = fluid.BOxmuO(pO).*bO;
    else
        muO = fluid.muO(pO);
    end
    mobO = krO./muO;
    dpO  = op.Grad(pO) - rhoOf.*gdz;
    upco = (double(dpO)<=0);
    
end