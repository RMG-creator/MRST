function [vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, dpW, ...
    extraOutput] = getFluxAndPropsWaterPolymer_BO(model, pO, sW, c, ...
    ads, krW, T, gdz)

f = model.fluid;
s = model.operators;

% Check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW') && ~isempty(sW)
    pcOW  = f.pcOW(sW);
end
pW = pO - pcOW;

% Multipliers due to polymer
mixpar   = f.mixPar;
cbar     = c./f.cmax;
a        = f.muWMult(f.cmax).^(1-mixpar);
b        = 1./(1-cbar+cbar./a);
muWMult  = f.muWMult(c); % viscosity multiplier from polymer mixing
muWMultT = b.*muWMult.^mixpar;

% Water props
bW      = f.bW(pO);
rhoW    = bW.*f.rhoWS;
rhoWf   = s.faceAvg(rhoW); % rhoW on face, average of neighboring cells
muW     = f.muW(pO);
muWeff  = muWMultT.*muW;
Rk      = 1 + ((f.rrf-1)./f.adsMax).*ads;
mobW    = krW./(muWeff.*Rk);
dpW     = s.Grad(pO-pcOW) - rhoWf.*gdz;
upcw    = (double(dpW)<=0); % water upstream-index
vW      = -s.faceUpstr(upcw, mobW).*T.*dpW;
if any(bW < 0)
    warning('Negative water compressibility present!')
end

% Polymer
mobP   = (mobW.*c)./(a + (1-a)*cbar);
vP     = - s.faceUpstr(upcw, mobP).*s.T.*dpW;

% Return extra output if requested
if model.extraPolymerOutput
    extraOutput.muWeff = muWeff;
    muPeff = muWeff.*(a + (1-a)*cbar);
    extraOutput.muPeff = muPeff;
else
    extraOutput = [];
end

end


