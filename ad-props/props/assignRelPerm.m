function f = assignRelPerm(f)

if ~isfield(f, 'krOG')      % two-phase water/oil
   if isfield(f, 'krWSurf')
      % Surfactant is used
      f.relPerm = @(sw, varargin)relPermWOSurf(sw, f, varargin{:});
   else
      f.relPerm = @(sw, varargin)relPermWO(sw, f, varargin{:});
   end
elseif ~isfield(f, 'krOW')  % two-phase oil/gas
    f.relPerm = @(sg, varargin)relPermOG(sg, f, varargin{:});
else                        % three-phase
    f.relPerm = @(sw, sg, varargin)relPermWOG(sw, sg, f, varargin{:});
end

end

function [krW, krO] = relPermWO(sw, f, varargin)
krW = f.krW(sw, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sw, varargin{:});
else
    krO = f.krOW(1-sw, varargin{:});
end
end


function [krW, krO, krWSurf, krOSurf] = relPermWOSurf(sw, f, varargin)
krW = f.krW(sw, varargin{:});
krWSurf = f.krWSurf(sw, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sw, varargin{:});
    krOSurf = f.krOSurf(1-sw, varargin{:});
else
    krO = f.krOW(1-sw, varargin{:});
    krOSurf = f.krOWSurf(1-sw, varargin{:});
end
end

function [krO, krG] = relPermOG(sg, f, varargin)
krG = f.krG(sg, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sg, varargin{:});
else
    krO = f.krOG(1-sg, varargin{:});
end
end

function [krW, krO, krG] = relPermWOG(sw, sg, f, varargin)
swcon = f.sWcon;
swcon = min(swcon, double(sw)-1e-5);

d  = (sg+sw-swcon);
ww = (sw-swcon)./d;
%krW = ww.*f.krW(sg+sw, varargin{:});
krW = f.krW(sw, varargin{:});

wg = 1-ww;
%krG = wg.*f.krG(sg+sw-swcon, varargin{:});
krG = f.krG(sg, varargin{:});

so = 1-sw-sg;
krow = f.krOW(so, varargin{:});
krog = f.krOG(so,  varargin{:});
krO  = wg.*krog + ww.*krow;
end
