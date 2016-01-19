%% Read case from file
% This example contains a simple $31\times31\times3$ fine grid containing
% two injectors in opposite corners and one producer in the middle of the
% domain. All wells are completed in the top layers of cells.
%
% The schedule being used contains first a period of injection with surfactant, followed by a water
% flooding phase without surfactant. Finally, the water rate is reduced for the final time steps.
%

try
   require add ad-core ad-blackoil ad-eor ad-props deckformat optimization
catch
   mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat optimization
end

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SURFACTANT.DATA');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);

gravity on

%% Set up simulation parameters
% We want a layer of oil on top of the reservoir and water on the bottom.
% To do this, we alter the initial state based on the logical height of
% each cell. The resulting oil concentration is then plotted.

ijk = gridLogicalIndices(G);

state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;

% Enforce s_w + s_o = 1;
state0.s(:,1) = 1 - state0.s(:,2);

% Add zero surfactant concentration to the state.
state0.c    = zeros(G.cells.num, 1);
state0.cmax = zeros(G.cells.num, 1);

clf
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
title('Oil concentration')
axis tight off
view(70, 30);
colorbar;



%% Set up systems.

modelSurfactant = OilWaterSurfactantModel(G, rock, fluid, 'inputdata', deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(G, modelSurfactant, rock, deck);


%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsSurfactant, statesSurfactant] = ...
   simulateScheduleAD(state0, modelSurfactant, schedule);

