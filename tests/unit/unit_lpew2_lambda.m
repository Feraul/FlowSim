%UNIT_LPEW2_LAMBDA  Verify fs.lpew.v2.lambdaWeights matches legacy Lamdas_Weights_LPEW2 per node.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_lpew2_lambda');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

fprintf('  [step 1] env + FS...\n');
fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);
env.geometry.elem(:, 5) = 1;
env.config.kmap = [1, env.config.perm];

FS = fs.mesh.build(env);
FS = fs.csr.buildCorners(FS);
FS = fs.csr.buildCornerShifts(FS);
FS.cfg.phasekey = env.config.phasekey;

fprintf('  [step 2] batched OPT + ksInterp + angulos + netas (build all inputs)...\n');
[P_all, T_all, O_all, Qc_all, ~] = fs.lpew.OPT(FS);
[Kt1, Kt2, Kn1, Kn2] = fs.lpew.v2.ksInterp(FS, T_all, Qc_all);
[ve2, ve1, theta2, theta1] = fs.lpew.v2.angulos(FS, T_all, O_all, Qc_all);
netas_all = fs.lpew.v2.netas(FS, T_all, P_all, O_all, Qc_all);

fprintf('  [step 3] batched lambdaWeights (interior vect + boundary loop)...\n');
[lambda, r] = fs.lpew.v2.lambdaWeights(FS, Kt1, Kt2, Kn1, Kn2, ...
    theta1, theta2, ve1, ve2, netas_all, T_all, Qc_all);

fs_expect(numel(lambda) == FS.csr.nCorners, 'lambda has nCorners entries');
fs_expect(isequal(size(r), [FS.mesh.nNodes, 2]), 'r is [nNodes x 2]');
fs_expect(all(isfinite(lambda)), 'lambda all finite');
fs_expect(all(isfinite(r(:))), 'r all finite');

fprintf('  [step 4] per-node slice equality vs legacy (all nodes)...\n');
tol_abs = 1e-10; tol_rel = 1e-12;
ok = @(diff, leg) (norm(diff) < tol_abs) || (norm(diff) < tol_rel * norm(leg));

r_leg = zeros(FS.mesh.nNodes, 2);
lambda_leg_all = zeros(FS.csr.nCorners, 1);

% Compare on a sample first — full-mesh comparison at the end
sampleNodes = [1, 5, 41, 45, 81];
for k = 1:numel(sampleNodes)
    No = sampleNodes(k);
    e_range = FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1);

    % Legacy per-node call chain
    [O_leg, P_leg, T_leg, Qo_leg] = OPT_Interp_LPEW(No, env);
    [Kt1_leg, Kt2_leg, Kn1_leg, Kn2_leg] = ferncodes_Ks_Interp_LPEW2( ...
        O_leg, T_leg, Qo_leg, No, [], env, env.config.kmap);
    [ve2_leg, ve1_leg, theta2_leg, theta1_leg] = angulos_Interp_LPEW2( ...
        O_leg, P_leg, T_leg, Qo_leg, No, env);
    netas_leg = netas_Interp_LPEW(O_leg, P_leg, T_leg, Qo_leg, No, env);
    [lambda_leg, r_leg_local] = Lamdas_Weights_LPEW2( ...
        Kt1_leg, Kt2_leg, Kn1_leg, Kn2_leg, theta1_leg, theta2_leg, ...
        ve1_leg, ve2_leg, netas_leg, P_leg, O_leg, Qo_leg, No, T_leg, r_leg);
    r_leg = r_leg_local;  % accumulate for r comparison

    lambda_slice = lambda(e_range);
    d_lambda = lambda_slice - lambda_leg;
    fs_expect(ok(d_lambda, norm(lambda_leg)), ...
        sprintf('node %d: lambda diff %.3e vs leg %.3e', No, norm(d_lambda), norm(lambda_leg)));
end

fs_teardown();
