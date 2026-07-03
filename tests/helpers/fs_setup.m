function fs_setup(testName)
%FS_SETUP  Standard test setup: addpath, seed rng, start timer, reset state.
%
%   fs_setup('smoke_env')  % logs test name, initialises fs_test_state
%
%   Adds FlowSim + its subfolders to path so tests written from anywhere
%   can resolve MetodoBase, ferncodes_*, etc.

    if nargin < 1, testName = 'unnamed'; end

    % Reset per-test state
    setappdata(0, 'fs_test_state', struct('pass', 0, 'fail', 0, 'first_fail', ''));
    setappdata(0, 'fs_test_name', testName);
    setappdata(0, 'fs_test_t0', tic);

    % Add FlowSim to path (pwd is expected to be the FlowSim repo root)
    for sub = {'', 'base', 'solvers', 'factories', 'simulacoes', 'benchmarks', ...
               'tests', 'tests/helpers', 'tests/smoke', 'tests/unit'}
        d = fullfile(pwd, sub{1});
        if isfolder(d), addpath(d); end
    end

    % Reproducibility
    rng(1234, 'twister');

    fprintf('\n╔══════════════════════════════════════════════════╗\n');
    fprintf('║  TEST: %-42s║\n', testName);
    fprintf('╚══════════════════════════════════════════════════╝\n');
end
