% UTILS
%
% Files
%   addFluxesFromSourcesAndBC       - Add in fluxes imposed by sources and face boundary conditions
%   assignValue                     - Assign values to ADI object by way of indices, without changing jacobians
%   bc2ADbc                         - Undocumented utility function
%   checkWellConvergence            - Compute convergence for wells.
%   CNV_MBConvergence               - Compute convergence based on total mass balance and maximum residual mass balance.
%   compressSchedule                - Compress schedule to take the longest possible timesteps while honoring controls
%   convertDeckScheduleToMRST       - Convert deck-type schedule to MRST style schedule
%   convertReportToSchedule         - Create a new schedule based on actual ministeps from a simulation report
%   convertReservoirFluxesToSurface - Compute surface fluxes from reservoir fluxes
%   double2ADI                      - Convert a double to ADI variable, using a sample ADI variable for dimensions
%   fastInterpTable                 - Undocumented utility function
%   fluidPlotPanelAD                - Create simple plot of fluid model for a given simulation model.
%   getBoundaryConditionFluxesAD    - Get boundary condition fluxes for a given set of values
%   getFaceTransmissibility         - Compute face transmissibilities, accounting for input-specific multipliers
%   getMultipliers                  - Get dynamic multiplier values for reservoir quantities
%   getPerforationToWellMapping     - Get map from global perforation number to global well index.
%   getSimulationTime               - Get the global time for a set of states produced by simulateScheduleAD
%   getSourceFluxesAD               - Short description
%   getWellOutput                   - Extract values from wellsols.
%   initWellSolAD                   - Set up well solution struct for a automatic differentiation model
%   makeScheduleConsistent          - Ensure that a schedule is consistent in terms of well counts/perforations
%   padRatesAndCompi                - Pad one/two/threephase values with zeros corresponding to missing phases.
%   pressureBCContrib               - Undocumented utility function
%   pressureBCContribADI            - Undocumented utility function
%   printConvergenceReport          - Undocumented utility function
%   recoverVars                     - Recover previously eliminated variables x at position n using solutions sol
%   reorderForILU                   - Undocumented utility function
%   ResultHandler                   - Class for storing and retrieving simulation results, either in memory or stored to disk
%   selectModelFromDeck             - Select simulation model from a ECLIPSE/FrontSim style input deck
%   setupOperatorsTPFA              - Set up helper structure for solvers based on automatic differentiation.
%   setWellSign                     - Ensure that wells have a defined sign. Will attempt to guess based on controls.
%   simpleSchedule                  - Make a schedule with varying timesteps and fixed wells/bc/src terms
%   terniaryWellPlot                - Undocumented utility function
%   uniqueStable                    - Support UNIQUE(A, 'stable') in all versions of MATLAB

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
