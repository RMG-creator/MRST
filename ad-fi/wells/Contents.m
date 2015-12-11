% WELLS
%
% Files
%   calcAlpha                      - Undocumented utility function
%   calcCompi                      - Undocumented utility function
%   checkLims                      - Undocumented utility function
%   computeWellContributions       - [eqs, cq_s, closedConns] = getWellContributions(...
%   computeWellHead                - compute pressure drop from ref depth to connections of corresponding well
%   eqsWellBO                      - Undocumented utility function
%   getActivePhases                - Undocumented utility function
%   getControlEquations            - Undocumented utility function
%   getWellContributions           - [eqs, cq_s, sol, closedConns] = getWellContributions(W, sol, pBH, q_s, p, rho_s, b, r, rMax, m, varargin)
%   getWellContributionsBO         - Undocumented utility function
%   reorderWellPerforationsByDepth - Undocumented utility function
%   solveLocalEqs                  - Solves the equations only for cells where the residual is large.
%   solveLocalWellEqs              - Undocumented utility function
%   updateConnDP                   - Explicit update of hydrostatic pressure difference between bottom hole
%   updateControls                 - Undocumented utility function
%   updateState                    - Undocumented utility function

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
