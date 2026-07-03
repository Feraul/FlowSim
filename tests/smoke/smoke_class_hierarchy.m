addpath(fullfile(pwd, 'tests', 'helpers'));  %% bootstrap so fs_* helpers resolve
%SMOKE_CLASS_HIERARCHY  Verify OOP class-hierarchy state (which load vs which fail).
%
%   Codifies the study finding: MetodoBase + SimulacaoBase are alive;
%   SolverBase + BenchmarkBase are missing (making SolverMPFAH, SolverNLFVPP,
%   Caso1 unloadable); 33 of 35 CasoNNN classes are missing.
%
%   This test EXPECTS the broken state. When the study eventually addresses
%   it (e.g. by creating SolverBase / BenchmarkBase / MetodoMPFAH / ...),
%   this test flips to expect the FIXED state — the failing assertions
%   become the "regression detector" for the fix.
fs_setup('smoke_class_hierarchy');

works = {'MetodoBase', 'MetodoMPFAD', 'MetodoTPFA', ...
         'SimulacaoBase', 'SimRichards', 'SimGroundwater', 'Caso439'};
for k = 1:numel(works)
    name = works{k};
    ok = false;
    try
        mc = meta.class.fromName(name);
        ok = ~isempty(mc);
    catch
        ok = false;
    end
    fs_expect(ok, sprintf('%s loads (should be alive)', name));
end

% Study finding: these classes should FAIL to load in current state.
% When study evolves, flip these expectations (e.g. add SolverBase.m → assert loads).
broken = {'SolverMPFAH', 'SolverNLFVPP', 'SolverBase', ...
          'BenchmarkBase', 'Caso1', 'Caso331', 'Caso437'};
for k = 1:numel(broken)
    name = broken{k};
    isEmpty = false; hasErr = false;
    try
        mc = meta.class.fromName(name);
        isEmpty = isempty(mc);
    catch
        hasErr = true;
    end
    fs_expect(isEmpty || hasErr, ...
        sprintf('%s does NOT load (study assertion: currently broken)', name));
end

fs_teardown();
