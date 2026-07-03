addpath(fullfile(pwd, 'tests', 'helpers'));  %% bootstrap so fs_* helpers resolve
%SMOKE_ENV  Sanity check: MATLAB version, cwd, path configuration.
fs_setup('smoke_env');

fs_expect(~verLessThan('matlab', '9.7'), ...
    'MATLAB >= R2019b (found R2024a expected)');

v = version;
fs_expect(contains(v, 'R2024a'), ...
    sprintf('MATLAB version is R2024a (got %s)', v));

fs_expect(isfolder(fullfile(pwd, 'tools')), ...
    'FlowSim tools/ folder present');

fs_expect(exist('tools/mrun', 'file') == 2, ...
    'tools/mrun is present');

fs_expect(isfolder(fullfile(pwd, 'meshes')) || any(~cellfun(@isempty, {dir('*.msh').name})), ...
    'mesh files reachable (either meshes/ subdir or root .msh)');

nRoot = numel(dir('*.m'));
fs_expect(nRoot > 100, sprintf('root .m file count is >100 (got %d — down from 285 pre-cleanup)', nRoot));

fs_expect(exist('preprocessormod', 'file') == 2, ...
    'preprocessormod is on path (canonical mesh preprocessor)');

fs_expect(exist('MetodoBase', 'class') == 8, ...
    'MetodoBase class resolves');

fs_teardown();
