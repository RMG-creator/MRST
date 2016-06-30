%% Buckley-Leverett problem with multiple substeps
% This example demonstrates how the sequential solvers perform on a one
% dimensional problem. In addition, we demonstrate how the sequential
% solvers can be configured to have a different time step selection for
% pressure and transport in order to reduce numerical diffusion from the
% temporal discretization.
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential

%% Set up model
% Construct 3D grid with 50 cells in the x-direction
G = cartGrid([50, 1, 1], [1000, 10, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', darcy*ones(G.cells.num, 1), ...
              'poro', .3*ones(G.cells.num, 1));

% Default oil-water fluid with unit values, quadratic relative permeability
% curves and a 2:1 viscosity ratio between oil and water.
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2],...
                           'mu', [1, 2]*centi*poise);

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 50*barsa, [0, 1]);
state0.wellSol = initWellSolAD([], model, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv = poreVolume(G, rock);
injRate = -sum(pv)/(500*day);
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);

%% Solve fully implicit base case
% We first solve the problem using the standard fully implicit
% discretization found in ad-blackoil in the regular manner.
dT = 20*day;
schedule = simpleSchedule(repmat(dT,1,25), 'bc', bc);
[~, states] = simulateScheduleAD(state0, model, schedule);
%% Get a sequential model
% We can set up a sequential model from the fully implicit model. The
% sequential model is a special model that solves each control step with a
% pressure and a transport model, each with their own linear and nonlinear
% solvers. The advantage of this abstraction is that it easy to create
% highly tailored solvers that exploit the properties of the underlying
% equations.
seqModel = getSequentialModelFromFI(model);
% We have two different models, each with their own properties and options,
% allowing for flexible configurations.
% Display the pressure model
disp(seqModel.pressureModel)
% Display the pressure model
disp(seqModel.transportModel);
%% Simulate sequential base case
% We can now simply pass the sequential model to simulateScheduleAD just
% like we did with the regular fully implicit model.
[~, statesSeq] = simulateScheduleAD(state0, seqModel, schedule);
%% Set up timestep selection based on target change quantities
% We will now compute a solution with refined time steps. As the time-steps
% becomes smaller, the solution becomes more accurate. In order to achieve
% increased accuracy without manually changing the timesteps, we can use
% a automatic timestep selector based on saturation change targets. We let
% the solver aim for a maximum saturation change of 1% in each cell during
% the timesteps to get very small steps.
stepSel = StateChangeTimeStepSelector(...
          'targetProps', 's',...            % Saturation as change target
          'targetChangeAbs', 0.01,...       % Target change of 0.01
          'firstRampupStepRelative', 0.01); % Initial rampup step is dt0/100

%% Simulate with refined timesteps
% The main issue with small timesteps is that they can be very time
% consuming to compute, especially as the full system has to be solved at
% every step. We therefore pass the step selector to the nonlinear solver
% corresponding to the transport subproblem in order to only take small
% steps in the transport solver itself, leaving the pressure only to be
% updated in the original control steps.
seqModel.transportNonLinearSolver.timeStepSelector = stepSel;
[~, statesFine, repFine] = simulateScheduleAD(state0, seqModel, schedule);
%% Plot and compare the three different solutions
% By plotting the three solutions simultanously, we see that the fully
% implicit and sequential implicit give the same solution for the same
% timesteps. The solution that is refined in time is quite different,
% however, demonstrating that the numerical truncation error due to the
% temporal discretization can significantly impact the simulation results.
for i = 1:numel(statesFine)
    figure(1); clf, hold on
    plot(states{i}.s(:, 1), '-')
    plot(statesSeq{i}.s(:, 1), '.')
    plot(statesFine{i}.s(:, 1), '--')
    ylim([0, 1])
    legend('Fully implicit', 'Sequential', 'Sequential \Delta s_{max} = 0.01');
    pause(0.1)
end