classdef OilWaterSurfactantModelAdsExpl < OilWaterSurfactantModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.

   properties
      explicitAdsComputation;
   end


   methods

      function model = OilWaterSurfactantModelAdsExpl(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantModel(G, rock, fluid, varargin{:});
         explicitAdsComputation = true;
      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
         [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                        'explicitAdsorption', true, varargin{:});
      end

      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

         [state, report] = updateAfterConvergence@OilWaterSurfactantModel(model, state0, state, ...
                                                           dt, drivingForces);

         state = updateExplicitAds(state0, state, model, dt);

      end


      function [fn, index] = getVariableField(model, name)
      % Get the index/name mapping for the model (such as where
      % pressure or water saturation is located in state)
         switch(lower(name))
           case {'ads'}
             index = 1;
             fn = 'ads';
           case {'adsmax'}
             index = 1;
             fn = 'adsmax';
           otherwise
             [fn, index] = getVariableField@OilWaterSurfactantModel(...
                model, name);
         end
      end

   end
end
