function [vW, vM, vN, bW, mobW, mobM, mobN, rhoW, pW, upcw]= ...
            getFluxAndPropsMEOR(model, pO, sW, m, n, krW, T, gdz)
fluid = model.fluid;
s = model.operators;

% Check for capillary pressure (p_cOW)
pcOW = 0;
if isfield(fluid, 'pcOW') && ~isempty(sW)
    pcOW = fluid.pcOW(sW);
end
pW = pO - pcOW;

% Will eventually edit to add viscosity changes and shit like that

bW = fluid.bW(pO);
rhoW = bW.*fluid.rhoWS;
% rhoW on face is the average of the neighboring cells
rhoWf = s.faceAvg(rhoW);
muW = fluid.muW(pO);
muWeff = muW;       %EDIT
mobW = krW./muWeff;
dpW = s.Grad(pO-pcOW) - rhoWf.*gdz;
% water upstream index
upcw = double(dpW)<=0;
vW = -s.faceUpstr(upcw,mobW).*T.*dpW;
if any(bW <0)
    warning('Negative water compressibility present')
end


% MEOR
mobM = mobW.*m;
mobN = mobW.*n;
vM = -s.faceUpstr(upcw, mobM).*s.T.*dpW;
vN = -s.faceUpstr(upcw, mobN).*s.T.*dpW;


