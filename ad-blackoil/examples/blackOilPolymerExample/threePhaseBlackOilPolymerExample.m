%%
% This example contains a simple 2D 20*5*1 grid containing one injector
% well on the bottom layer and one production well on the top layer.
%
% The schedule being used contains first a period of injection without
% polymer, then water flooding period with polymer, followed by a pure
% water flooding phase without polymer.
%

clear;

mrstModule add ad-core ad-props ad-blackoil ad-fi deckformat

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, '2D_THREEPHASE_POLY_HETER.DATA');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

%% set the initial conditsions with with InitEclipseState
% for this deck, equil is used to initialize the initial state, while the
% results of the initialization need to be verified against Eclipse later.

gravity on;

fluid          = initDeckADIFluid(deck);
state0         = initEclipseState(G, deck, initEclipseFluid(deck));
state0.c       = zeros([G.cells.num, 1]);
state0.cmax    = zeros([G.cells.num, 1]);
modelBOPolymer = ThreePhaseBlackOilPolymerModel(G, rock, fluid, 'inputdata', deck);

modelBOPolymer.disgas = 1;
modelBOPolymer.vapoil = 1;

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(G, modelBOPolymer, rock, deck);

%% Set the non-linear solver
modelBOPolymer.useCNVConvergence = true;

if ~isempty(mrstPath('agmg'))
    mrstModule add agmg
    pSolver = AGMGSolverAD();
    linsolve = CPRSolverAD('ellipticSolver', pSolver);
    nonlinearsolver = NonLinearSolver('LinearSolver', linsolve);
    nonlinearsolver.useRelaxation = true
else
    nonlinearsolver = NonLinearSolver();
    nonlinearsolver.useRelaxation = true;
end

%% plotting the initial saturation
h=figure(1);
set(h, 'Position', [100, 100, 900, 600]);
clf;

subplot(2,2,1);
plotGrid(G, 'facec', 'none');
axis tight off;
view(26,18);
W = schedule.control(1).W;
plotWell(G,W);

subplot(2,2,2);
plotCellData(G, state0.s(:,1));
plotGrid(G, 'facec', 'none')
title('Initial Water saturation')
axis tight off
view(26,18);
colorbar;

subplot(2,2,3);
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
title('Initial Oil saturation')
axis tight off
view(26,18);
colorbar;

subplot(2,2,4);
plotCellData(G, state0.s(:,3));
plotGrid(G, 'facec', 'none')
title('Initial Gas saturation')
axis tight off
view(26,18);
colorbar;

pause(0.5);

%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.


[wellSolsPolymer, statesPolymer] = ...
   simulateScheduleAD(state0, modelBOPolymer, schedule, 'NonLinearSolver', nonlinearsolver);


%% plotting the wells data
h = figure(3);
set(h, 'Position', [100, 50, 900, 1300]);
clf;

T     = convertTo(cumsum(schedule.step.val), day);

[qWs, qOs, qGs, bhp] = wellSolToVector(wellSolsPolymer);

subplot(5,2,1);
plot(T, bhp(:,1)/barsa);
title('BHPs for injection wells');
ylabel('bhp (bar)');
xlabel('time (day)');

subplot(5,2,2);
plot(T, bhp(:,2)/barsa);
title('BHPs for production wells');
ylabel('bhp (bar)');
xlabel('time (day)');

subplot(5,2,3);
plot(T, qWs(:,1)*day);
title('water injection rate');
ylabel('WWIR (m^3/day)');
xlabel('time (day)');

subplot(5,2,4);
plot(T, qWs(:,2)*day);
title('Water Production Rate');
ylabel('WWPR (m^3/day)');
xlabel('time (day)');

subplot(5,2,5);
plot(T, qOs(:,1)*day);
title('oil injection rate');
ylabel('WOIR (m^3/day)');
xlabel('time (day)');

subplot(5,2,6);
plot(T, qOs(:,2)*day);
title('oil production rate');
ylabel('WOPR (m^3/day)');
xlabel('time (day)');

subplot(5,2,7);
plot(T, qGs(:,1)*day);
title('gas injection rate');
ylabel('WGIR (m^3/day)');
xlabel('time (day)');

subplot(5,2,8);
plot(T, qGs(:,2)*day);
title('gas production rate');
ylabel('WGPR (m^3/day)');
xlabel('time (day)');


extract1 = @(wsol, fld) [ wsol.(fld) ];
extract2 = @( c ) vertcat(c{:});
extract  = @(fld) extract2(cellfun(@(wsol) extract1(wsol, fld), ...
                           wellSolsPolymer, 'UniformOutput', false));

sgn =        extract('sign');
qWPoly = sgn .* extract('qWPoly');

subplot(5,2,9);
plot(T, qWPoly(:,1)*day);
title('polymer injection rate');
ylabel('WCIR (kg/day)');
xlabel('time (day)');

subplot(5,2,10);
plot(T, qWPoly(:,2)*day);
title('polymer production rate');
ylabel('WCPR (kg/day)');
xlabel('time (day)');

pause(0.5);


%% plotting the animated saturation evolution
numStep = size(statesPolymer, 1);
for iStep = 1:10:numStep

    h = figure(2);
    set(h, 'Position', [700, 100, 900, 600]);
    clf;


    stateStep = statesPolymer{iStep,1};

    subplot(2,2,1);
    plotCellData(G, stateStep.s(:,1));
    plotGrid(G, 'facec', 'none')
    title('Water saturation')
    axis tight off
    view(26,18);
    colorbar;


    subplot(2,2,2);
    plotCellData(G, stateStep.s(:,2));
    plotGrid(G, 'facec', 'none')
    title('Oil saturation')
    axis tight off
    view(26,18);
    colorbar;

    subplot(2,2,3);
    plotCellData(G, stateStep.s(:,3));
    plotGrid(G, 'facec', 'none')
    title('Gas saturation')
    axis tight off
    view(26,18);
    colorbar;

    subplot(2,2,4);
    plotCellData(G, stateStep.c);
    plotGrid(G, 'facec', 'none')
    title('Polymer concentration')
    axis tight off
    view(26,18);
    colorbar;

    pause(0.5);

end
%%
save resMRSTPolymer wellSolsPolymer statesPolymer schedule;
fprintf('The simulation has been finished! \n');
