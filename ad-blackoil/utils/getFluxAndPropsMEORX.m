function [vW, vM_0, vM_1, vPol, vN, bW, mobW, mobM_0, mobM_1, mobPol, mobN, rhoW, pW, upcw]= ...
            getFluxAndPropsMEORX(model, pO, sW, m_0, m_1, pol, n, krW, T, gdz)
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
mobM_0 = mobW.*m_0;
mobM_1 = mobW.*m_1;
mobPol = mobW.*pol;
mobN = mobW.*n;
vM_0 = -s.faceUpstr(upcw, mobM_0).*s.T.*dpW;
vM_1 = -s.faceUpstr(upcw, mobM_1).*s.T.*dpW;
vPol = -s.faceUpstr(upcw, mobPol).*s.T.*dpW;
vN = -s.faceUpstr(upcw, mobN).*s.T.*dpW;


