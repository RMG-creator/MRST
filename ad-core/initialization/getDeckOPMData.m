function [deck, output] = getDeckOPMData(name, caseName)
    opm = mrstPath('opm-data');
    opm_tests = mrstPath('opm-tests');
    fn = fullfile(opm, name, [caseName, '.DATA']);
    deck = readEclipseDeck(fn);
    
    units = {'METRIC', 'FIELD', 'SI', 'LAB', 'PVT_M'};
    
    unit = 'SI';
    for i = 1:numel(units)
        if isfield(deck.RUNSPEC, units{i})
            unit = units{i};
            break
        end
    end
    
    deck = convertDeckUnits(deck);
    fldr = fullfile(opm_tests, name);
    output = struct('opm', [], 'ecl', []);
    try
        output.opm = getOutputOPM(deck, fldr, caseName, unit);
    catch e
        fprintf('Error in parsing of OPM output. Reported error:\n %s\n',e.message);
    end
    try
        output.eclipse = getOutputEclipse(deck, fldr, caseName, unit);
    catch e
        fprintf('Error in parsing of Eclipse output. Reported error:\n %s\n',e.message);
    end
    
end


function output = getOutputOPM(deck, fldr, caseName, unit)
    G = initEclipseGrid(deck, 'SplitDisconnected', false);
    smry_prefix = fullfile(fldr, 'opm-simulation-reference', 'flow', caseName);
    if ~exist([smry_prefix, '.UNSMRY'], 'file')
        smry_prefix = fullfile(fldr, 'opm-simulation-reference', 'flow_legacy', caseName);
    end

    [ws, wstime] = convertSummaryToWellSols(smry_prefix, unit);
    states = convertRestartToStates(smry_prefix, G, 'use_opm', true, 'consistentWellSols', false);

    output = struct('wellSols', {ws}, 'states', {states}, 'G', G, 'location', smry_prefix, 'time', wstime);
end

function output = getOutputEclipse(deck, fldr, caseName, unit)
    G = initEclipseGrid(deck);
    smry_prefix = fullfile(fldr, 'eclipse-simulation', caseName);

    ws = convertSummaryToWellSols(smry_prefix, unit);
    states = convertRestartToStates(smry_prefix, G, 'use_opm', false, 'consistentWellSols', false);

    output = struct('wellSols', {ws}, 'states', {states}, 'G', G);
end