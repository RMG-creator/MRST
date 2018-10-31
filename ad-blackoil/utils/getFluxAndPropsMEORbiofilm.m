function [vW, vMicro, vMeta, vN, bW, mobW, mobm, mobMeta, mobN, rhoW, pW, upcw]= ...
            getFluxAndPropsMEORbiofilm(model, pO, sW, m, meta, n, bio, krW, T, gdz)
fluid = model.fluid;
s = model.operators;
% hack for avoiding complex values
inx = meta<0;
meta(inx) = 0;

% Check for capillary pressure (p_cOW)
pcOW = 0;
if isfield(fluid, 'pcOW') && ~isempty(sW)
    pcOW = fluid.pcOW(sW);
end
pW = pO - pcOW;

change = 1.*ones(length(meta),1); %0 for power law 1 for parabolic
if model.biopoly
    % Threshold for metabolite effect/hack to avoid complex values
    inx = meta<1e-16;
    meta_eff = meta;
    meta_eff(inx) = 0;
    % power law from Lacerda
    %change = 1.4019.*meta_eff.^.1653;
    % Parabolic law from Bartelds
    change = ((5.*meta_eff).^2 + 5.*meta_eff + 1);
    inx = change<1e-16;
    change(inx) = 0;
end

bW = fluid.bW(pO);
rhoW = bW.*fluid.rhoWS;
% rhoW on face is the average of the neighboring cells
rhoWf = s.faceAvg(rhoW);
muW = fluid.muW(pO);
%muWeff = muW + change.*1e-3; % For power law      
muWeff = muW.*change;          % For parabolic law
%disp(max(double(muWeff)))
mobW = krW./muWeff;
dpW = s.Grad(pO-pcOW) - rhoWf.*gdz;
% water upstream index
upcw = double(dpW)<=0;
vW = -s.faceUpstr(upcw,mobW).*T.*dpW;
if any(bW <0)
    warning('Negative water compressibility present')
end


% MEOR with biofilm
mobm = mobW.*(m);
mobMeta = mobW.*meta;
mobN = mobW.*n;
vMicro = -s.faceUpstr(upcw, mobm).*T.*dpW;
vMeta = -s.faceUpstr(upcw, mobMeta).*T.*dpW;
vN = -s.faceUpstr(upcw, mobN).*T.*dpW;


