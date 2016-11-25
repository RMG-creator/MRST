classdef TwoPhaseOilGasModel < ThreePhaseBlackOilModel
   % Two phase oil/gas system
   properties
   end

   methods
      function model = TwoPhaseOilGasModel(G, rock, fluid, varargin)
         model = ...
            model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

         % This is the model parameters for oil/gas
         model.oil   = true;
         model.gas   = true;
         model.water = false;

         % Blackoil -> use CNV style convergence
         model.useCNVConvergence = true;

         model.saturationVarNames = { 'so', 'sg' };
         model.wellVarNames       = { 'qOs', 'qGs', 'bhp' };

         model = merge_options(model, varargin{:});
      end

      %--------------------------------------------------------------------

      function [problem, state] = ...
            getEquations(model, state0, state, dt, drivingForces, varargin)

         [problem, state] = ...
            equationsBlackOil(state0, state, model, dt, ...
                              drivingForces, varargin{:});
      end
   end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
