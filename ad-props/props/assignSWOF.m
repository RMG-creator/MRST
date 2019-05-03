function f = assignSWOF(f, swof, reg)
    [f.krW, f.krOW, f.pcOW, f.krPts.w, f.krPts.ow] = getFunctions(swof, reg);
end

function [krW, krOW, pcOW, pts_w, pts_ow] = getFunctions(SWOF, reg)
    pts_w = zeros(reg.sat, 4);
    pts_ow = zeros(reg.sat, 4);
    
    [krW, krOW, pcOW] = deal(cell(1, reg.sat));
    
    for i = 1:reg.sat
        [pts_w(i, :), pts_ow(i, :)] = getPoints(SWOF{i});
        
        swof = extendTab(SWOF{i});
        SW = swof(:, 1);
        krW{i}  = @(sw) interpTable(SW, swof(:, 2), sw);
        krOW{i} = @(so) interpTable(SW, swof(:, 3), 1-so);
        pcOW{i} = @(sw) interpTable(SW, swof(:, 4), sw);
    end
end

function [pts, pts_o] = getPoints(swof)
    pts = zeros(1, 4);
    % Connate water saturation
    pts(1) = swof(1, 1);
    % Last immobile water saturation
    ii = find(swof(:,2)==0, 1, 'last');
    pts(2) = swof(ii,1);
    % Last point
    pts(3) = swof(end,1);
    % Maximum relperm
    pts(4) = swof(end,2);
    
    % Get OW-scaling
    pts_o = zeros(1, 4);
    pts_o(3) = 1;
    ii = find(swof(:,3) == 0, 1, 'first');
    pts_o(2) = 1 - swof(ii,1);
    % Maximum oil relperm
    pts_o(4) = swof(1,3);
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

