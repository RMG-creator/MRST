function [krW_eff , krO_eff , krG_eff , muW_eff , muO_eff , muG_eff , rhoW_eff, rhoO_eff, rhoG_eff, bW_eff  , bO_eff  , bG_eff  , bW0_eff , bO0_eff , bG0_eff , pvMult, transMult, mobMult, pvMult0, T] =  getDynamicQuantitiesOilWaterSolvent(model, p0, p, sW, sO, sG, sO0, sG0)

    fluid = model.fluid;
    op    = model.operators;
    
    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    
    
    %% Transmissibility
    T = op.T.*transMult;
    
    %% Relative permeabilites
    
    % Residual saturations
    sOres = fluid.sOres;
    sGres = fluid.sGres;
    
%     sOn = max(sO - sOres,0);
%     sGn = max(sG - sGres,0);
    
    sOn = sO;
    sGn = sG;

    sNn = sOn + sGn;
    
%     sOn0 = max(sO0 - sOres,0);
%     sGn0 = max(sG0 - sGres,0);
%     
    sOn0 = sO0;
    sGn0 = sG0;
    
    sNn0 = sOn0 + sGn0;
    
    % Effective relative permeabilites
    krN = @(s) fluid.krO(s);
    krW_eff = fluid.krW(sW);

    krO_eff = sOn./sNn.*krN(sO + sG);
    krG_eff = sGn./sNn.*krN(sO + sG);
   
    % Multiply by mobility multiplyer
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
        
    %% Viscosities
    
%     % Unmixed phase viscosities (water not affected by solvent)
%     muO = fluid.muO(p);
%     muG = fluid.muG(p);
%     
%     % Mixing parameter
%     omega = fluid.mixPar;
%     
%     % Mixed viscosity
%     a = 1/4;
%     muM = muO.*muG./(sGn./sNn.*muO.^a + sOn./sNn.*muG.^a).^4;
%     
%     % Effective viscosities
%     muW_eff = fluid.muW(p);
%     muO_eff = muO.^(1-omega).*muM.^omega;
%     muG_eff = muG.^(1-omega).*muM.^omega;

    [muW, muO, muG, muW_eff , muO_eff , muG_eff ] = calculateViscosities(fluid, p, sOn, sGn, sNn);
    [~  , muO0  , ~  , muW_eff0, muO_eff0, muG_eff0] = calculateViscosities(fluid, p0, sOn0, sGn0, sNn0);
    
    %% Densities
    
%     % Unmixed phase viscosities (water not affected by solvent)
%     rhoO = fluid.bO(p).*fluid.rhoOS;
%     rhoG = fluid.bG(p).*fluid.rhoGS;
%     
%     
%     if any(abs(muO - muG) > eps)
%         % Expressions are valid if muO ~= muG
%         
%         % Effective satuiration fractions
%         r = muO./muG;
%         sR_Oeff = (r.^a - (muO./muO_eff).^a)./(r.^a - 1);
%         sR_Geff = (r.^a - (muO./muG_eff).^a)./(r.^a - 1);
% 
%         % Effective densities
%         rhoW_eff = fluid.bW(p).*fluid.rhoWS;
%         rhoO_eff = rhoO.*sR_Oeff + rhoG.*(1 - sR_Oeff);
%         rhoG_eff = rhoO.*sR_Geff + rhoG.*(1 - sR_Geff);
%     
%     else
%         % Expressions are sinuglar if muO == muG
%         
%         rhoM = rhoO.*(sO./sN) + rhoG.*(sG./sN);
%         rhoO_eff = (1-omega).*rhoO + omega.*rhoM;
%         rhoG_eff = (1-omega).*rhoG + omega.*rhoM;
%         
%     end

[rhoW_eff, rhoO_eff, rhoG_eff] = calculateDensities(fluid, p, muO_eff, muG_eff, sOn, sGn, sNn); 
[rhoW_eff0, rhoO_eff0, rhoG_eff0] = calculateDensities(fluid, p0, muO_eff0, muG_eff0, sOn, sGn, sNn); 
        
    %% Formation volume factors
    
    bW_eff = fluid.bW(p);
%     bO = fluid.bO(p);
    bO_eff = rhoO_eff./fluid.rhoOS;
    bG_eff = rhoG_eff./fluid.rhoGS;
%     bG = fluid.bG(p);
        
    bW0_eff = fluid.bW(p0);
%     bO0 = fluid.bO(p0);
    bO0_eff = rhoO_eff0./fluid.rhoOS;
%     bG0 = fluid.bG(p0);
    bG0_eff = rhoG_eff0./fluid.rhoGS;
        
end 

function [muW, muO, muG, muW_eff, muO_eff, muG_eff] = calculateViscosities(fluid, p, sOn, sGn, sNn)

 % Unmixed phase viscosities (water not affected by solvent)
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    
    % Mixing parameter
    omega = fluid.mixPar;
    
    % Mixed viscosity
    a = 1/4;
    muM = muO.*muG./(sGn./sNn.*muO.^a + sOn./sNn.*muG.^a).^4;
    
    % Effective viscosities
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muM.^omega;
    muG_eff = muG.^(1-omega).*muM.^omega;
    
end

function [rhoW_eff, rhoO_eff, rhoG_eff] = calculateDensities(fluid, p, muO_eff, muG_eff, sO, sG, sN)

    % Unmixed phase viscosities (water not affected by solvent)
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    
    muO = fluid.muO(p);
    muG = fluid.muG(p);
        
    omega = fluid.mixPar;
    
    a = 1/4;
    % Expressions are valid if muO ~= muG

    % Effective satuiration fractions
    r = muO./muG;
    sR_Oeff = (r.^a - (muO./muO_eff).^a)./(r.^a - 1);
    sR_Geff = (r.^a - (muO./muG_eff).^a)./(r.^a - 1);

    % Effective densities
    rhoW_eff = fluid.bW(p).*fluid.rhoWS;
    
    
    % Expressions are sinuglar if muO == muG

    rhoM = rhoO.*(sO./sN) + rhoG.*(sG./sN);
    
    tol = 1e-10;
    eq = abs(muO - muG) < tol;
    
    rhoO_eff = rhoO.*sR_Oeff + rhoG.*(1 - sR_Oeff).*(~eq) + (1-omega).*rhoO + omega.*rhoM.*eq;
    rhoG_eff = rhoO.*sR_Geff + rhoG.*(1 - sR_Geff).*(~eq) + (1-omega).*rhoG + omega.*rhoM.*eq;
    
%     rhoO_eff = ;
%     rhoG_eff = (1-omega).*rhoG + omega.*rhoM;
    
    
    
end