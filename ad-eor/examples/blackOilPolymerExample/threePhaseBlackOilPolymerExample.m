%% 2D Three-Phase Polymer Injection Case
% This example contains a simple 4000 m × 200 m × 125 m reservoir described
% on 20 × 1 × 5 uniform Cartesian grid. One injection well is located in
% the bottom two layers and one production well is located in the top two
% layers. Hydrostatic equilibration is used for initialization.

% The polymer injection schedule follows a typical polymer waterflooding
% strategy. The flooding process begins with primary waterflooding for 1260
% days, followed by a polymer slug injected over 1700 days, and then
% switching back to water injection. The injection well is under rate
% control with target rate 1000 m3 /day and upper limit of 450 bar on the
% bottom-hole pressure (bhp), whereas the production well is under pressure
% control with target bottom-home pressure 260 bar.

mrstModule add ad-core ad-blackoil ad-eor ad-fi ad-props deckformat

%% Set up model and initial conditions
% The data required for the example
% If the data does not exist locally, download it automatically
fname = { '2D_THREEPHASE_POLY_HETER.DATA', 'POLY.inc' };
files = fullfile(getDatasetPath('BlackoilPolymer2D', 'download', true), fname);

% check to make sure the files are complete
e = cellfun(@(pth) exist(pth, 'file') == 2, files);

if ~all(e),
    pl = ''; if sum(e) ~= 1, pl = 's'; end
    msg = sprintf('Missing data file%s\n', pl);
    msg = [msg, sprintf('  * %s\n', fname{~e})];
    error('Dataset:Incomplete', msg);
end

% Parsing the data file
deck = readEclipseDeck(files{1});
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

gravity on;
fluid          = initDeckADIFluid(deck);
modelBOPolymer = ThreePhaseBlackOilPolymerModel(G, rock, fluid, 'inputdata', deck);

% Set initial equilibrium and no polymer present
state0         = initEclipseState(G, deck, initEclipseFluid(deck));
state0.c       = zeros([G.cells.num, 1]);
state0.cmax    = zeros([G.cells.num, 1]);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(modelBOPolymer, deck);

%% Select nonlinear and linear solvers
% Use the AGMG algebraic multigrid solver if this is present, and if not,
% use the default MATLAB solver.
modelBOPolymer.useCNVConvergence = true;

% Choose the linear solver for pressure equation
if ~isempty(mrstPath('agmg'))
    mrstModule add agmg
    pSolver = AGMGSolverAD();
else
    pSolver = BackslashSolverAD();
end
linsolve = CPRSolverAD('ellipticSolver', pSolver);
nonlinearsolver = getNonLinearSolver(modelBOPolymer, 'DynamicTimesteps', false, 'LinearSolver', linsolve);
nonlinearsolver.useRelaxation = true;


%% Plotting the initial saturations
h=figure(1); clf
set(h, 'Position', [100, 100, 900, 600]);

subplot(2,2,1);
plotGrid(G, 'FaceColor', 'none'); plotWell(G,schedule.control(1).W);
axis tight off; view(26,18);

phases = {'water','oil','gas'};
for i=1:3
    subplot(2,2,i+1);
    plotCellData(G, state0.s(:,i), 'EdgeColor','k');
    title(['Initial ', phases{i}, ' saturation'])
    axis tight off; view(26,18); colorbar;
end

%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsPolymer, statesPolymer] = ...
    simulateScheduleAD(state0, modelBOPolymer, schedule, ...
                       'NonLinearSolver', nonlinearsolver);

%% plotting the well data
figure(2)
set(gcf,'Position', [100, 50, 900, 1300]);
T = convertTo(cumsum(schedule.step.val), day);
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSolsPolymer);

subplot(5,2,1); title('BHPs for injection wells');
plot(T, convertTo(bhp(:,1), barsa));
ylabel('bhp (bar)'); xlabel('time (day)');

subplot(5,2,2); title('BHPs for production wells');
plot(T, convertTo(bhp(:,2), barsa));
ylabel('bhp (bar)'); xlabel('time (day)');

subplot(5,2,3); title('water injection rate');
plot(T, convertTo(qWs(:,1), meter^3/day));
ylabel('WWIR (m^3/day)'); xlabel('time (day)');

subplot(5,2,4); title('Water Production Rate');
plot(T, convertTo(qWs(:,2), meter^3/day));
ylabel('WWPR (m^3/day)'); xlabel('time (day)');

subplot(5,2,5); title('oil injection rate');
plot(T, convertTo(qOs(:,1), meter^3/day));
ylabel('WOIR (m^3/day)'); xlabel('time (day)');

subplot(5,2,6); title('oil production rate');
plot(T, convertTo(qOs(:,2), meter^3/day));
ylabel('WOPR (m^3/day)'); xlabel('time (day)');

subplot(5,2,7); title('gas injection rate');
plot(T, convertTo(qGs(:,1), meter^3/day));
ylabel('WGIR (m^3/day)'); xlabel('time (day)');

subplot(5,2,8); title('gas production rate');
plot(T, convertTo(qGs(:,2), meter^3/day));
ylabel('WGPR (m^3/day)'); xlabel('time (day)');

qWPoly = getWellOutput(wellSolsPolymer, 'qWPoly');
sign = getWellOutput(wellSolsPolymer, 'sign');
qWPoly = sign .* qWPoly;

subplot(5,2,9); title('polymer injection rate');
plot(T, convertTo(qWPoly(:,1), kilogram/day));
ylabel('WCIR (kg/day)'); xlabel('time (day)');

subplot(5,2,10); title('polymer production rate');
plot(T, convertTo(qWPoly(:,2), kilogram/day));
ylabel('WCPR (kg/day)'); xlabel('time (day)');

%% plotting the animated saturation evolution
numStep = size(statesPolymer, 1);
figure(2); set(gcf,'Position', [700, 100, 900, 600]);
for iStep = 1:10:numStep
    stateStep = statesPolymer{iStep,1};

    for i=1:3,
        subplot(2,2,i), cla
        plotCellData(G, stateStep.s(:,i),'EdgeColor','k');
        title([phases{i} ' saturation']);
        axis tight off, view(26,18); colorbar; caxis([0 1]);
    end

    subplot(2,2,4); cla
    plotCellData(G, stateStep.c,'EdgeColor','k');
    title('Polymer concentration')
    axis tight off; view(26,18); colorbar; caxis([0 1]);

    drawnow; pause(0.5);

end
%%
% save resMRSTPolymer wellSolsPolymer statesPolymer schedule;
fprintf('The simulation has been finished! \n');
