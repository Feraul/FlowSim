function [env, weight, s] = preLPEW2(env, parms)
%FS.LPEW.V2.PRELPEW2  Vectorized LPEW2 driver (drop-in replacement for ferncodes_Pre_LPEW_2_vect).
%
%   [env, weight, s] = fs.lpew.v2.preLPEW2(env, parms)
%
%   Runs the full LPEW2 pipeline as batched ops:
%     1. Build FS + CSR + shifts
%     2. Batched OPT → P, T, O, Qcorner
%     3. Batched ksInterp → Kt1/Kt2/Kn1/Kn2 (all corners)
%     4. Batched angulos + netas → theta1/theta2, ve1/ve2, netas (all corners)
%     5. Batched lambdaWeights → lambda (interior vect + boundary loop), r
%     6. Segmented normalisation:  weight = lambda / sum(lambda per node)
%     7. Neumann source s via env.benchmark.calcularTermoNeumannVet
%
%   Same output layout as ferncodes_Pre_LPEW_2_vect: weight is a 1×nCorners
%   row-vector aligned with esurn1 flat order. s is [nNodes x 1].

    % ── Ensure the kmap the caller expected ────────────────────────
    if env.config.numcase > 400 && isfield(parms, 'auxperm') && ~isempty(parms.auxperm)
        env.config.kmap = parms.auxperm;
    elseif ~isfield(env.config, 'kmap') || isempty(env.config.kmap)
        env.config.kmap = env.config.perm;
    end

    % ── Build FS scaffolding ───────────────────────────────────────
    FS = fs.mesh.build(env);
    FS = fs.csr.buildCorners(FS);
    FS = fs.csr.buildCornerShifts(FS);
    FS.cfg.phasekey = env.config.phasekey;

    % ── Batched kernels ────────────────────────────────────────────
    [P_all, T_all, O_all, Qc_all, ~] = fs.lpew.OPT(FS); %#ok<ASGLU>
    [Kt1, Kt2, Kn1, Kn2] = fs.lpew.v2.ksInterp(FS, T_all, Qc_all);
    [ve2, ve1, theta2, theta1] = fs.lpew.v2.angulos(FS, T_all, O_all, Qc_all);
    netas_all = fs.lpew.v2.netas(FS, T_all, P_all, O_all, Qc_all);
    [lambda, r] = fs.lpew.v2.lambdaWeights(FS, Kt1, Kt2, Kn1, Kn2, ...
        theta1, theta2, ve1, ve2, netas_all, T_all, Qc_all);

    % ── Segmented normalisation: weight = lambda / sum(lambda per node) ─
    sum_lambda = accumarray(FS.csr.cornerNode, lambda, [FS.mesh.nNodes, 1]);
    % Guard against 0/0 (isolated nodes with no lambda contribution — shouldn't
    % happen on a well-formed mesh but be defensive)
    denom = sum_lambda(FS.csr.cornerNode);
    denom(denom == 0) = 1;
    weight_col = lambda ./ denom;

    % Row-vector layout matching legacy Pre_LPEW_2_vect output
    weight = weight_col.';

    % ── Neumann source term (benchmark-provided) ───────────────────
    N = env.premethod.MPFAD.N;
    s = env.benchmark.calcularTermoNeumannVet(r, sum_lambda, N, env);

    % ── Persist in env.premethod ───────────────────────────────────
    env.premethod.MPFAD.weight = weight;
    env.premethod.MPFAD.s      = s;
end
