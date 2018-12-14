%% sensitivitiesModel2D - analyse sensitivity capabilities 

mrstModule add ad-core ad-blackoil ad-props optimization spe10 mrst-gui

% Setup model -> grid, rock, schedule, fluid etc
setupModel2D

sc = schedule;
sc.step.val = sc.step.val(1 : 3);
sc.step.control = sc.step.control(1 : 3);
schedule = sc;

%% Reset fluid to include scaling:
% $s_w -> \frac{s_w-swcr}{swu-swcr}$
% $s_o -> \frac{s_o-sowcr}{1-swl-sowcr}$

% fluid = initSimpleScaledADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                                 % 'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                                 % 'n',     [2, 2, 0], ...
                                 % 'swl',   0.10*ones(G.cells.num,1), ...
                                 % 'swcr',  0.15*ones(G.cells.num,1), ...
                                 % 'sowcr', 0.12*ones(G.cells.num,1), ...
                                 % 'swu',   0.90*ones(G.cells.num,1));

fluid = initSimpleScaledADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                                 'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                                 'n',     [2, 2, 0]);

                       
% Create model-object of class TwoPhaseOilWaterModel
model_ref  = TwoPhaseOilWaterModel(G, rock, fluid);                       

% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [.15, .85]);

% Set up a perturbed model with different pv and perm:
rock1 = rock;
rock1.perm = rock.perm*1.1;
model = TwoPhaseOilWaterModel(G, rock1, fluid);                       
model.operators.pv = model_ref.operators.pv.*0.8;

% run ref model
[ws_ref, states_ref, r_ref] = simulateScheduleAD(state0, model_ref, schedule);
% run model
[ws, states, r] = simulateScheduleAD(state0, model, schedule);

% plot well solutions for the two models
% plotWellSols({ws_ref, ws}, {r_ref.ReservoirTime, r.ReservoirTime}, ...
            % 'datasetnames', {'reference', 'perturbed'})

%% setup misfit-function and run adjoint to get parameter sensitivities
% setup weights for matching function, empty weight uses default (will 
% produce function value of ~O(1) for 100% misfit). Only match rates in this example: 
weighting =  {'WaterRateWeight',     [], ...
              'OilRateWeight',       [] , ...
              'BHPWeight',           0};
   
% compute misfit function value (first each summand corresonding to each time-step)
misfitVals = matchObservedOW(G, ws, schedule, ws_ref, weighting{:});
% sum values to obtiain scalar objective 
misfitVal = sum(vertcat(misfitVals{:}));
fprintf('Current misfit value: %6.4e\n', misfitVal)

% setup (per time step) mismatch function handle for passing on to adjoint sim
objh = @(tstep) matchObservedOW(G, ws, schedule, ws_ref, 'computePartials', ...
                                true, 'tstep', tstep, weighting{:});
% objhm = @(tstep) matchObservedOW(G, m.ws, schedule, m.ws_ref, 'computePartials', ...
                                % true, 'tstep', tstep, weighting{:});

params      = {'transmissibility', 'conntrans'};
paramTypes  = {'multiplier', 'multiplier'};
sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , params, ...
                                     'ParameterTypes', paramTypes);


figure,
subplot(1,2,1), plotCellData(G,  cellAverage(G, sens.transmissibility), 'EdgeColor', 'none'), colorbar,title('Average trans multiplier sensitivity');
subplot(1,2,2), plot(sens.conntrans, '--o', 'MarkerSize', 14); title('Well connection trans multiplier sensitivity')
a = gca; a.XTick = 1:4;  a.XTickLabel = {W.name}; a.XLim; a.XLim = [.5 4.5];


