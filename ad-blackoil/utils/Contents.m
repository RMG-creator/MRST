% UTILS
%
% Files
%   assignWellValuesFromControl        - Utility function to assign well values that are the resut of the controls
%   calculateHydrocarbonsFromStatusBO  - define sG, rs and rv in terms of x
%   computeFlashBlackOil               - Update state based on vapoil/disgas flash
%   equationsBlackOil                  - Generate linearized problem for a volatile 3Ph system (wet-gas, live-oil).
%   equationsOilWater                  - Get linearized problem for oil/water system with black oil-style
%   equationsOilWaterPolymer           - Get linearized problem for oil/water/polymer system with black oil-style
%   equationsThreePhaseBlackOilPolymer - Get linearized problem for oil/water/gas/polymer system with black
%   equationsWater                     - Get linearized problem for single phase water system with black oil-style
%   equationsWaterThermal              - Get linearized problem for single phase water system with black oil-style
%   getbG_BO                           -
%   getbO_BO                           -
%   getCapillaryPressureBO             -
%   getCellStatusVO                    - Status should be passed on from updateStateVO (to be sure definition is
%   getFluxAndPropsGas_BO              -
%   getFluxAndPropsOil_BO              -
%   getFluxAndPropsWater_BO            -
%   getFluxAndPropsWaterPolymer_BO     -
%   updateStateBlackOilGeneric         - Generic update function for blackoil-like models

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
