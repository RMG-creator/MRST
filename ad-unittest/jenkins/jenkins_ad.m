cd('../../../')
ls
run('mrst-core/startup.m');

% Set up module directories
names = {'autodiff', ...
         'internal', ...
         'multiscale', ...
         'visualization', ...
         'model-io', ...
         'solvers', ...
         'thirdparty-modules'};
     
names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);
mrstPath('addroot', names{:});

names = {'test-datasets'};
names = cellfun(@(x) fullfile(ROOTDIR, '..', x), names, ...
                    'UniformOutput', false);
mrstPath('register', 'co2lab', fullfile(ROOTDIR, '..', 'mrst-co2lab', 'co2lab'));
mrstPath('register', 'agmg', fullfile('..', 'hnil-agmg'));
mrstPath('addroot', names{:});

% Ensure Matlab BGL is available
bgl_dir = fullfile(ROOTDIR(), 'utils', '3rdparty','matlab_bgl');
if ~exist(fullfile(bgl_dir, 'matlab_bgl'), 'dir')
    run(fullfile(bgl_dir,'downloadMBGL'))
end


% Remove previous tests
jfolder = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'XML');
if exist(jfolder, 'dir')
    rmdir(jfolder, 's')
end

tap_folder = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'JUnit');
if exist(tap_folder, 'dir')
    rmdir(tap_folder, 's')
end

mrstModule add ad-unittest ad-core ad-blackoil
runTestsAD('runIntegration',    true, ...
            'runUnit',          true, ...
            'runExamples',      true, ...
            'writeToDisk',      true, ...
            'writeXML',         true);
% Disable matlab
exit();