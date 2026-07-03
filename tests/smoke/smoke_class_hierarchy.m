addpath(fullfile(pwd, 'tests', 'helpers'));  %% bootstrap so fs_* helpers resolve
%SMOKE_CLASS_HIERARCHY  Verify OOP class-hierarchy is unified under MetodoBase.
%
%   Post-PR-A2 state:
%     - MetodoBase + subclasses MPFAD, TPFA, MPFAH, NLFVPP, MPFAQL: all load
%     - SimulacaoBase + subclasses SimRichards, SimGroundwater, Caso439: all load
%     - SolverBase, BenchmarkBase, SolverMPFAH, SolverNLFVPP, Caso1: gone/dead
%     - Missing CasoNNN (331, 437, ...): still absent (out of PR-A2 scope)
fs_setup('smoke_class_hierarchy');

works = {'MetodoBase', 'MetodoMPFAD', 'MetodoTPFA', ...
         'MetodoMPFAH', 'MetodoNLFVPP', 'MetodoMPFAQL', ...
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
    fs_expect(ok, sprintf('%s loads (unified MetodoBase hierarchy)', name));
end

% Post-PR-A2: these must NOT load — they were deleted / never should have existed.
gone = {'SolverBase', 'SolverMPFAH', 'SolverNLFVPP', 'BenchmarkBase', 'Caso1'};
for k = 1:numel(gone)
    name = gone{k};
    isEmpty = false; hasErr = false;
    try
        mc = meta.class.fromName(name);
        isEmpty = isempty(mc);
    catch
        hasErr = true;
    end
    fs_expect(isEmpty || hasErr, ...
        sprintf('%s does NOT load (deleted in PR-A2)', name));
end

% Missing CasoNNN classes — still absent (Phase F triage).
stillMissing = {'Caso331', 'Caso437'};
for k = 1:numel(stillMissing)
    name = stillMissing{k};
    try
        mc = meta.class.fromName(name);
        isEmpty = isempty(mc);
    catch
        isEmpty = false;
    end
    fs_expect(isEmpty, sprintf('%s still absent (Phase F scope)', name));
end

fs_teardown();

