function flowsim_init(varargin)
%FLOWSIM_INIT  Add all FlowSim source folders to the MATLAB path.
%
%   flowsim_init                  % default: run silently, verify version
%   flowsim_init('verbose', true) % print every folder added
%   flowsim_init('legacy', false) % skip legacy/ (fully-vectorized runs only)
%   flowsim_init('reset',  true)  % rmpath the tree first, then re-add
%
%   Ordering:
%     1. +fs/**                (new vectorized tree — resolves FIRST)
%     2. benchmarks/, factories/, base/, solvers/   (transitional)
%     3. legacy/**             (old procedural code — resolves LAST so
%                               vectorized twins shadow it during migration)
%
%   Run this ONCE per MATLAB session. Symmetric cleanup: flowsim_deinit.
%
%   Author: FlowSim vectorization project (AXON code-dev flowsim-vectorize)
%   Date:   2026-07-03

    % ── argument parsing ─────────────────────────────────────────────────
    p = inputParser;
    addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'legacy',  true,  @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'reset',   false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    opt = p.Results;

    % ── locate repo root (this file MUST live at repo root) ──────────────
    rootDir = fileparts(mfilename('fullpath'));

    % ── MATLAB version check (>= R2019b for +package support + arguments block) ─
    v = version('-release');
    year = str2double(v(1:4));
    if year < 2019
        error('flowsim_init:MatlabTooOld', ...
              'FlowSim requires MATLAB R2019b or newer (found R%s).', v);
    end

    % ── optional reset: remove previously-added FlowSim paths ────────────
    if opt.reset
        flowsim_deinit('verbose', opt.verbose);
    end

    % ── collect paths in the correct precedence order ────────────────────
    paths = {};

    % 1. new vectorized tree (+fs/**) — resolves first
    fsRoot = fullfile(rootDir, '+fs');
    if isfolder(fsRoot)
        paths{end+1} = genpath(fsRoot); %#ok<AGROW>
    else
        warning('flowsim_init:MissingFs', ...
                '+fs/ not found at %s — vectorized tree unavailable.', fsRoot);
    end

    % 2. transitional OOP layer (already env-based, no globals)
    for sub = {'benchmarks', 'factories', 'base', 'solvers', 'simulacoes'}
        d = fullfile(rootDir, sub{1});
        if isfolder(d), paths{end+1} = d; end %#ok<AGROW>
    end

    % 3. legacy procedural code — resolves last (shadowable by +fs twins)
    if opt.legacy
        legacyRoot = fullfile(rootDir, 'legacy');
        if isfolder(legacyRoot)
            paths{end+1} = genpath(legacyRoot); %#ok<AGROW>
        else
            % Migration not yet done — legacy files still live at root.
            % Include repo root itself, but *after* everything above.
            paths{end+1} = rootDir; %#ok<AGROW>
        end
    end

    % ── add to MATLAB path in reverse order so the FIRST entry has highest
    %    precedence after addpath processes them (addpath prepends by default)
    %    We want +fs/ prepended LAST → highest priority.
    for k = numel(paths):-1:1
        addpath(paths{k}, '-end');   % add legacy at end
    end
    % Now re-add +fs at the front to guarantee shadowing:
    if isfolder(fsRoot)
        addpath(genpath(fsRoot), '-begin');
    end

    % ── verbose summary ──────────────────────────────────────────────────
    if opt.verbose
        fprintf('╔════════════════════════════════════════════════════════════╗\n');
        fprintf('║  FlowSim path initialized                                  ║\n');
        fprintf('║  Root:     %-48s║\n', rootDir);
        fprintf('║  Verbose:  true                                            ║\n');
        fprintf('║  Legacy:   %-48s║\n', ternary(opt.legacy, 'included', 'skipped'));
        fprintf('╚════════════════════════════════════════════════════════════╝\n');
        for k = 1:numel(paths)
            fprintf('  + %s\n', paths{k});
        end
    else
        fprintf('FlowSim ready (MATLAB R%s, legacy=%d). Run: flowsim_run\n', ...
                v, opt.legacy);
    end

    % ── stash init state so deinit can undo cleanly ──────────────────────
    setappdata(0, 'flowsim_init_paths', paths);
    setappdata(0, 'flowsim_init_root',  rootDir);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
