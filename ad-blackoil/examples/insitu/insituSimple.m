%% Gravity segregation using two phase AD solvers
mrstVerbose on
mrstModule add ad-core ad-blackoil ad-props ad-fi mrst-gui
[schedule, model, state0] = getBenchmarkAD('spe9');

% lim = 10;
% schedule.step.val = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);

%% Simulate the problem
close all
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, 'plotReservoir', false);

mrstModule add agmg
linsolve = CPRSolverAD('ellipticSolver', AGMGSolverAD());
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
    'Verbose', true, 'afterStepFn', fn, 'linearsolver', linsolve, 'outputministeps', true);