% startup.m — auto-invoked by MATLAB at session start when cwd is FlowSim root.
% Now delegates to the canonical flowsim_init.m (PR-A1) which handles all
% paths (base/simulacoes/benchmarks/solvers/factories + +fs/ + legacy/) with
% correct precedence.
function startup()
    base = fileparts(mfilename('fullpath'));
    initFile = fullfile(base, 'flowsim_init.m');

    if exist(initFile, 'file') == 2
        % Prefer the canonical initializer — handles +fs/ + legacy/ + more
        run(initFile);
        return;
    end

    % Fallback: legacy startup behaviour (kept for pre-PR-A1 checkouts)
    addpath(fullfile(base, 'base'));
    addpath(fullfile(base, 'simulacoes'));
    addpath(fullfile(base, 'benchmarks'));
    addpath(fullfile(base, 'solvers'));
    addpath(fullfile(base, 'factories'));
    fprintf('Paths configurados com sucesso.\n');
end