function [FS, capture] = capture_baseline(varargin)
%CAPTURE_BASELINE  Run FlowSim end-to-end on a fixture mesh + capture the numerical fingerprint.
%
%   [FS, capture] = capture_baseline('mesh', 'M8.msh', 'method', 'tpfa', 'numcase', 439)
%
%   Runs the reduced main.m pipeline:
%     1. fs_test_env — isolate Start.dat into a WSL-writable temp
%     2. preprocessormod(1) — build env from patched Start.dat
%     3. createBenchmark / createMetodo / createSimulacao
%     4. sim.preprocessar + PLUG_kfunction + ferncodes_calflag + preprocessmethod
%     5. env.metodo.montarSistema  → M, I  (STOP HERE for baseline — no solver yet)
%
%   Saves capture struct to tests/golden/<mesh>-num<N>-<method>.mat with fields:
%     capture.mesh, capture.method, capture.numcase
%     capture.nNodes, capture.nElems, capture.nBFaces, capture.nIFaces
%     capture.M_size, capture.M_nnz, capture.M_frobnorm
%     capture.I_norm
%     capture.premethod_keys, capture.premethod_frobs  (per-field snapshot)
%     capture.ts, capture.matlab_version
%
%   Also returns the FS struct (from fs.mesh.build) for downstream diffing.
%
%   Called by tests/unit/unit_baseline_reproduces.m — which captures + diffs.

    % ── argument parsing ─────────────────────────────────────────────
    p = inputParser;
    addParameter(p, 'mesh',     'M8.msh',  @ischar);
    addParameter(p, 'method',   '',        @ischar);   % override pmethod ('tpfa'|'mpfad')
    addParameter(p, 'numcase',  [],        @isnumeric);
    addParameter(p, 'outfile',  '',        @ischar);
    parse(p, varargin{:});

    % ── locate FlowSim root ──────────────────────────────────────────
    fsRoot = getappdata(0, 'fs_test_root');
    if isempty(fsRoot), fsRoot = pwd; end

    % ── 1. isolated env + preprocessormod ────────────────────────────
    fs_test_env('mesh', p.Results.mesh);
    env = preprocessormod(1);

    % Optional override of pmethod / numcase (post-Start.dat)
    if ~isempty(p.Results.method),  env.config.pmethod  = p.Results.method;  end
    if ~isempty(p.Results.numcase), env.config.numcase  = p.Results.numcase; end

    % ── 2. instantiate the three objects ─────────────────────────────
    env.benchmark = createBenchmark(env.config.numcase);
    env.metodo    = createMetodo(env.config.pmethod);
    sim           = createSimulacao(env.config.phasekey);

    % ── 3. pre-processing pipeline ───────────────────────────────────
    parms        = env.benchmark.initParms();
    [env, parms] = sim.preprocessar(env, parms);
    [env, parms] = PLUG_kfunction(env, parms, 0);
    [env]        = ferncodes_calflag(env, parms, 0);
    [env, parms] = preprocessmethod(env, parms);

    % ── 4. one system assembly (no solver) ───────────────────────────
    dt = 1;
    try
        [M, I] = env.metodo.montarSistema(env, parms, dt);
        assemblyOK = true;
    catch err
        warning('capture_baseline:AssemblyFailed', ...
                'montarSistema failed: %s — capturing pre-assembly state only', err.message);
        M = []; I = [];
        assemblyOK = false;
    end

    % ── 5. build FS struct via the new adapter ───────────────────────
    FS = fs.mesh.build(env);
    FS = fs.csr.buildCorners(FS);

    % ── 6. compose the capture struct ────────────────────────────────
    capture = struct();
    capture.ts             = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
    capture.matlab_version = version;
    capture.fs_root        = fsRoot;
    capture.mesh           = p.Results.mesh;
    capture.method         = env.config.pmethod;
    capture.numcase        = env.config.numcase;
    capture.phasekey       = env.config.phasekey;

    % Mesh sizes
    capture.nNodes  = size(env.geometry.coord,  1);
    capture.nElems  = size(env.geometry.elem,   1);
    capture.nBFaces = size(env.geometry.bedge,  1);
    capture.nIFaces = size(env.geometry.inedge, 1);
    capture.nCorners = FS.csr.nCorners;

    % Assembly fingerprint
    capture.assemblyOK = assemblyOK;
    if assemblyOK
        capture.M_size     = size(M);
        capture.M_nnz      = nnz(M);
        capture.M_frobnorm = norm(M, 'fro');
        capture.M_max      = full(max(abs(M(:))));
        capture.I_norm     = norm(I);
        capture.I_max      = max(abs(I));
    end

    % premethod snapshot (structural + Frobenius per numeric field)
    premethodKey = upper(env.config.pmethod);
    if strcmp(premethodKey, 'MPFAD'), premethodKey = 'MPFAD'; end
    if strcmp(premethodKey, 'TPFA'),  premethodKey = 'TPFA';  end
    if isfield(env.premethod, premethodKey)
        pm = env.premethod.(premethodKey);
        fns = fieldnames(pm);
        capture.premethod_keys  = fns;
        capture.premethod_frobs = cell(numel(fns), 1);
        for k = 1:numel(fns)
            v = pm.(fns{k});
            if isnumeric(v) && ~isempty(v)
                capture.premethod_frobs{k} = norm(v(:));
            else
                capture.premethod_frobs{k} = NaN;
            end
        end
    end

    % ── 7. save ──────────────────────────────────────────────────────
    outfile = p.Results.outfile;
    if isempty(outfile)
        stem = sprintf('%s-num%g-%s', strrep(p.Results.mesh, '.msh', ''), ...
                       capture.numcase, capture.method);
        outfile = fullfile(fsRoot, 'tests', 'golden', [stem '.mat']);
    end
    outDir = fileparts(outfile);
    if ~isfolder(outDir), mkdir(outDir); end
    save(outfile, 'capture', '-v7');

    fprintf('captured → %s\n', outfile);
    fprintf('  mesh=%s method=%s numcase=%g phasekey=%g\n', ...
            capture.mesh, capture.method, capture.numcase, capture.phasekey);
    fprintf('  nNodes=%d nElems=%d nBFaces=%d nIFaces=%d nCorners=%d\n', ...
            capture.nNodes, capture.nElems, capture.nBFaces, capture.nIFaces, capture.nCorners);
    if assemblyOK
        fprintf('  M: %dx%d nnz=%d frob=%.6e  |  I norm=%.6e\n', ...
                capture.M_size(1), capture.M_size(2), capture.M_nnz, ...
                capture.M_frobnorm, capture.I_norm);
    else
        fprintf('  (assembly failed — mesh+premethod snapshot only)\n');
    end
end
