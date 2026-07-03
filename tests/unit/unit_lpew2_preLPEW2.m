%UNIT_LPEW2_PRELPEW2  End-to-end: fs.lpew.v2.preLPEW2 vs legacy ferncodes_Pre_LPEW_2_vect.
%
%   Runs the full FlowSim pipeline (preprocessormod → createBenchmark →
%   preprocessmethod), then calls BOTH the legacy and new preLPEW2
%   drivers on the same env and asserts weight+s bit-tolerance match.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_lpew2_preLPEW2');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

fprintf('  [step 1] full pipeline setup (as capture_baseline)...\n');
fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);
env.config.numcase = 439;   % Caso439 wired
env.config.pmethod = 'mpfad';

env.benchmark = createBenchmark(env.config.numcase);
env.metodo    = createMetodo(env.config.pmethod);
sim           = createSimulacao(env.config.phasekey);
parms         = env.benchmark.initParms();
[env, parms]  = sim.preprocessar(env, parms);
[env, parms]  = PLUG_kfunction(env, parms, 0);
[env]         = ferncodes_calflag(env, parms, 0);
[env, parms]  = preprocessmethod(env, parms);
% preprocessmethod already ran ferncodes_Pre_LPEW_2_vect and stored the result
weight_leg = env.premethod.MPFAD.weight;
s_leg      = env.premethod.MPFAD.s;

fprintf('  [step 2] call fs.lpew.v2.preLPEW2 on the same env...\n');
[env2, weight_new, s_new] = fs.lpew.v2.preLPEW2(env, parms); %#ok<ASGLU>

fs_expect(isequal(size(weight_leg), size(weight_new)), ...
    sprintf('weight size matches (legacy %s vs new %s)', mat2str(size(weight_leg)), mat2str(size(weight_new))));
fs_expect(isequal(size(s_leg), size(s_new)), 's size matches');

% Full-mesh tolerance
absTol = 1e-10; relTol = 1e-12;
ok = @(diff, leg) (norm(diff) < absTol) || (norm(diff) < relTol * norm(leg));

d_w = weight_new - weight_leg;
d_s = s_new - s_leg;

fs_expect(ok(d_w, norm(weight_leg)), ...
    sprintf('weight full-mesh diff %.3e vs leg norm %.3e', norm(d_w), norm(weight_leg)));
fs_expect(ok(d_s, norm(s_leg)), ...
    sprintf('s full-mesh diff %.3e vs leg norm %.3e', norm(d_s), norm(s_leg)));

% Per-element weight check (element ordering must match)
fs_expect(all(isfinite(weight_new)), 'weight finite');
fs_expect(all(isfinite(s_new)),      's finite');
fs_expect(all(weight_new >= -1e-10), 'weight all >= 0 (partition of unity)');

% Partition of unity check per node
FS = fs.mesh.build(env);
FS = fs.csr.buildCorners(FS);
sums = accumarray(FS.csr.cornerNode, weight_new(:), [FS.mesh.nNodes, 1]);
% Interior nodes should sum to 1; boundary nodes may not
FS = fs.csr.buildCornerShifts(FS);
interior_sums = sums(FS.csr.nodeIsInterior);
if ~isempty(interior_sums)
    max_dev = max(abs(interior_sums - 1));
    fs_expect(max_dev < 1e-10, sprintf('interior partition of unity holds (max dev %.3e)', max_dev));
end

fs_teardown();
