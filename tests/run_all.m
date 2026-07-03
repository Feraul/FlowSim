%RUN_ALL  Run all smoke + unit tests, aggregate results, exit 0/1.
%
%   Usage:  tools/mrun -c $(pwd) tests/run_all.m
%
%   Doesn't use exit() per-subtest — collects results into a matrix and
%   only exits at the very end. If any test fails, exit code is 1.

fprintf('\n=== FlowSim test harness ===\n');
fprintf('  pwd: %s\n', pwd);
fprintf('  ts:  %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

testFiles = [...
    dir(fullfile(pwd, 'tests', 'smoke', 'smoke_*.m')); ...
    dir(fullfile(pwd, 'tests', 'unit',  'unit_*.m'))];

results = struct('name', {}, 'pass', {}, 'total', {}, 'elapsed', {}, 'status', {});
totalFail = 0;

% run each test in a subprocess so exit() doesn't kill the runner
tmpDir = fullfile(tempdir, 'fs_tests');
if ~isfolder(tmpDir), mkdir(tmpDir); end

for k = 1:numel(testFiles)
    tf = testFiles(k);
    stem = regexprep(tf.name, '\.m$', '');
    outFile = fullfile(tmpDir, [stem '.txt']);

    % Fork: use system to run each test in a fresh matlab.exe -batch.
    % This isolates exit() calls and prevents state bleed.
    cmd = sprintf(['matlab.exe -batch "addpath(''%s''); addpath(''%s''); ' ...
                   'addpath(''%s''); addpath(''%s''); ' ...
                   'diary(''%s''); try; run(''%s''); diary off; catch e; ' ...
                   'fprintf(2, ''UNCAUGHT: %%s\\n'', e.message); diary off; exit(1); end"'], ...
        strrep(pwd, '\', '\\'), ...
        strrep(fullfile(pwd, 'base'), '\', '\\'), ...
        strrep(fullfile(pwd, 'tests', 'helpers'), '\', '\\'), ...
        strrep(fullfile(pwd, fileparts(tf.folder)), '\', '\\'), ...
        strrep(outFile, '\', '\\'), ...
        strrep(fullfile(tf.folder, tf.name), '\', '\\'));
    [status, ~] = system(cmd);

    r.name    = stem;
    r.status  = (status == 0);
    if ~r.status, totalFail = totalFail + 1; end
    r.elapsed = 0;
    if exist(outFile, 'file'), fprintf('%s\n', fileread(outFile)); end
    results(end+1) = r;
end

fprintf('\n=== summary ===\n');
for k = 1:numel(results)
    r = results(k);
    fprintf('  %s  %s\n', ternchar(r.status, '[OK]', '[FAIL]'), r.name);
end
fprintf('\n%d tests, %d failed\n', numel(results), totalFail);
exit(double(totalFail > 0));

function s = ternchar(cond, a, b)
    if cond, s = a; else, s = b; end
end
