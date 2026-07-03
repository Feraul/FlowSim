%UNIT_ASSEMBLY_MPFAD  Verify fs.assembly.mpfad.build matches legacy globalmatrix_MPFAD.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_assembly_mpfad');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

fprintf('  [step 1] full pipeline through preprocessmethod...\n');
fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);
env.config.numcase = 439;
env.config.pmethod = 'mpfad';
env.benchmark = createBenchmark(env.config.numcase);
env.metodo    = createMetodo(env.config.pmethod);
sim           = createSimulacao(env.config.phasekey);
parms         = env.benchmark.initParms();
[env, parms]  = sim.preprocessar(env, parms);
[env, parms]  = PLUG_kfunction(env, parms, 0);
[env]         = ferncodes_calflag(env, parms, 0);
[env, parms]  = preprocessmethod(env, parms);

fprintf('  [step 2] legacy assembly...\n');
[M_leg, I_leg, ~] = ferncodes_globalmatrix_MPFAD(env, parms);
fprintf('     legacy M: %dx%d nnz=%d frob=%.6e  I norm=%.6e\n', ...
    size(M_leg,1), size(M_leg,2), nnz(M_leg), norm(M_leg,'fro'), norm(I_leg));

fprintf('  [step 3] vectorized assembly...\n');
[M_new, I_new, ~] = fs.assembly.mpfad.build(env, parms);
fprintf('     new    M: %dx%d nnz=%d frob=%.6e  I norm=%.6e\n', ...
    size(M_new,1), size(M_new,2), nnz(M_new), norm(M_new,'fro'), norm(I_new));

fs_expect(isequal(size(M_leg), size(M_new)), 'M size matches');
fs_expect(isequal(size(I_leg), size(I_new)), 'I size matches');
fs_expect(nnz(M_leg) == nnz(M_new), sprintf('M nnz matches (%d)', nnz(M_leg)));

% Frobenius diff on M
frob_diff = norm(M_leg - M_new, 'fro');
frob_leg  = norm(M_leg, 'fro');
rel_frob  = frob_diff / max(1e-30, frob_leg);
fs_expect(rel_frob < 1e-12, sprintf('M Frobenius rel diff %.3e < 1e-12', rel_frob));

% L2 diff on I
l2_diff = norm(I_leg - I_new);
l2_leg  = norm(I_leg);
if l2_leg > 1e-15
    rel_i = l2_diff / l2_leg;
    fs_expect(rel_i < 1e-10, sprintf('I L2 rel diff %.3e < 1e-10', rel_i));
else
    fs_expect(l2_diff < 1e-10, sprintf('I abs diff %.3e < 1e-10 (leg is zero)', l2_diff));
end

% Also verify vs golden
try
    gold = load(fullfile(fsRoot, 'tests', 'golden', 'M8-num439-mpfad.mat'));
    gold = gold.capture;
    fs_expect(nnz(M_new) == gold.M_nnz, sprintf('M nnz matches golden (%d)', gold.M_nnz));
    frob_gold_diff = abs(norm(M_new, 'fro') - gold.M_frobnorm) / max(1e-30, gold.M_frobnorm);
    fs_expect(frob_gold_diff < 1e-12, sprintf('M Frobenius matches golden (rel %.3e)', frob_gold_diff));
catch err
    fprintf('  (skipping golden compare: %s)\n', err.message);
end

fs_teardown();
