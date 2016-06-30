%% Read case from file
% This example contains a simple $31\times31\times3$ fine grid containing two
% injectors in opposite corners and one producer in the middle of the
% domain. All wells are completed in the top layers of cells.
%
% The schedule being used contains first a period of injection with surfactant,
% followed by a water flooding phase without surfactant. Finally, the water rate
% is reduced for the final time steps.
%

try
    require ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui
catch
    mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui
end

current_dir = fileparts(mfilename('fullpath'));
simul_case = '1D'; % '1D' or '2D'

switch simul_case
  case '1D'
    fn = fullfile(current_dir, 'SURFACTANT1D.DATA');
    gravity off
  case '2D'
    fn = fullfile(current_dir, 'SURFACTANT.DATA');
    gravity on
  otherwise
    error('simul_case not recognized.');
end

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);


G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);



%% Set up simulation parameters
% We want a layer of oil on top of the reservoir and water on the bottom.
% To do this, we alter the initial state based on the logical height of
% each cell. The resulting oil concentration is then plotted.

switch simul_case
  case '1D'

    nc = G.cells.num;
    state0 = initResSol(G, 300*barsa, [ .2, .8]);

    % Add zero surfactant concentration to the state.
    state0.c    = zeros(G.cells.num, 1);
    state0.cmax = state0.c;

  case '2D'

    ijk = gridLogicalIndices(G);

    state0 = initResSol(G, 300*barsa, [ .9, .1]);
    state0.s(ijk{3} == 1, 2) = .9;
    state0.s(ijk{3} == 2, 2) = .8;

    % Enforce s_w + s_o = 1;
    state0.s(:,1) = 1 - state0.s(:,2);

    % Add zero surfactant concentration to the state.
    state0.c    = zeros(G.cells.num, 1);
    state0.cmax = state0.c;

    clf
    plotCellData(G, state0.s(:,2));
    plotGrid(G, 'facec', 'none')
    title('Oil concentration')
    axis tight off
    view(70, 30);
    colorbar;

  otherwise
    error('simul_case not recognized.')
end


%% Set up systems.

modelSurfactant = FullyImplicitOilWaterSurfactantModel(G, rock, fluid, ...
                                                  'inputdata', deck, ...
                                                  'extraStateOutput', true);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(modelSurfactant, deck);


%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

state0.ads = computeEffAds(state0.c, 0, modelSurfactant.fluid);
state0.adsmax = state0.ads;

resulthandler = ResultHandler('dataDirectory', pwd, 'dataFolder', 'cache', 'cleardir', true);
[wellSolsSurfactant, statesSurfactant] = simulateScheduleAD(state0, modelSurfactant, ...
                                                  schedule, 'OutputHandler', ...
                                                  resulthandler);

switch simul_case
  case '1D'
    plotToolbar(G, statesSurfactant, 'startplayback', true, 'plot1d', true)
  case '2D'
    plotToolbar(G, statesSurfactant, 'startplayback', true)
  otherwise
    error('simul_case not recognized.');
end