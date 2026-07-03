%UNIT_LPEW2_KSINTERP  Batched fs.lpew.v2.ksInterp vs legacy per-node ferncodes_Ks_Interp_LPEW2.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_lpew2_ksInterp');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

fprintf('  [step 1] building env from M8.msh via fs_test_env + preprocessormod...\n');
fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);
fprintf('  [step 1 done] nNodes=%d nElems=%d\n', size(env.geometry.coord,1), size(env.geometry.elem,1));

% Normalize kmap + elem material ids for isolated test:
% preprocessormod leaves env.config.perm = [K11 K12 K21 K22] (1x4, single material)
% but elem(:,5) can be 1..nElems (physical tags from gmsh). Legacy production
% code overwrites elem(:,5) via PLUG_kfunction. For this isolation test we
% force a single-material world: elem(:,5)=1 everywhere + kmap = [id K11..K22].
env.geometry.elem(:, 5) = 1;
env.config.kmap = [1, env.config.perm];   % [1 x 5]  cols 2..5 = tensor

fprintf('  [step 2] building FS + CSR + shifts...\n');
FS = fs.mesh.build(env);
FS = fs.csr.buildCorners(FS);
FS = fs.csr.buildCornerShifts(FS);
FS.cfg.phasekey = env.config.phasekey;

fprintf('  [step 3] batched OPT (all corners)...\n');
[P_all, T_all, O_all, Qc_all, ~] = fs.lpew.OPT(FS); %#ok<ASGLU>

fprintf('  [step 4] batched ksInterp (all corners)...\n');
[Kt1, Kt2, Kn1, Kn2] = fs.lpew.v2.ksInterp(FS, T_all, Qc_all);

fprintf('  [step 5] shape assertions...\n');
fs_expect(isequal(size(Kt1), [FS.csr.nCorners, 2]), 'Kt1 is [nCorners x 2]');
fs_expect(isequal(size(Kn1), [FS.csr.nCorners, 2]), 'Kn1 is [nCorners x 2]');
fs_expect(size(Kt2, 1) == FS.csr.nCorners && size(Kt2, 2) == 1, 'Kt2 is [nCorners x 1]');
fs_expect(size(Kn2, 1) == FS.csr.nCorners && size(Kn2, 2) == 1, 'Kn2 is [nCorners x 1]');
fs_expect(all(isfinite(Kt1(:))) && all(isfinite(Kn1(:))), 'Kt1/Kn1 all finite');
fs_expect(all(isfinite(Kt2(:))) && all(isfinite(Kn2(:))), 'Kt2/Kn2 all finite');

fprintf('  [step 6] per-node slice equality vs legacy (5 sample nodes)...\n');
sampleNodes = [1, 5, 41, 45, 81];
tol = 1e-12;
Sw = [];  % single-phase — Sw unused by legacy path

for k = 1:numel(sampleNodes)
    No = sampleNodes(k);
    fprintf('    node %d: computing legacy per-node + comparing slice...\n', No);

    [O_leg, P_leg, T_leg, Qo_leg] = OPT_Interp_LPEW(No, env); %#ok<ASGLU>
    [Kt1_leg, Kt2_leg, Kn1_leg, Kn2_leg] = ferncodes_Ks_Interp_LPEW2( ...
        O_leg, T_leg, Qo_leg, No, Sw, env, env.config.kmap);

    e_range = FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1);
    Kt1_slice = Kt1(e_range, :);
    Kt2_slice = Kt2(e_range);
    Kn1_slice = Kn1(e_range, :);
    Kn2_slice = Kn2(e_range);

    fs_expect(isequal(size(Kt1_slice), size(Kt1_leg)), sprintf('node %d: Kt1 size matches', No));
    fs_expect(isequal(size(Kn1_slice), size(Kn1_leg)), sprintf('node %d: Kn1 size matches', No));

    % Use absolute-OR-relative tolerance — for isotropic K, Kt1/Kt2 should be
    % exactly 0. Legacy has fp noise ~1e-27, ours is exact 0. Relative diff
    % is spuriously 1.0 in that case. Absolute floor (1e-10) catches this.
    absTol = 1e-10;
    relTol = 1e-12;
    ok = @(diff, leg) (norm(diff) < absTol) || (norm(diff) < relTol * norm(leg));

    d_Kt1 = Kt1_slice - Kt1_leg;  n_Kt1 = norm(Kt1_leg, 'fro');
    d_Kt2 = Kt2_slice - Kt2_leg;  n_Kt2 = norm(Kt2_leg);
    d_Kn1 = Kn1_slice - Kn1_leg;  n_Kn1 = norm(Kn1_leg, 'fro');
    d_Kn2 = Kn2_slice - Kn2_leg;  n_Kn2 = norm(Kn2_leg);

    fs_expect(ok(d_Kt1(:), n_Kt1), sprintf('node %d: Kt1 diff %.3e vs leg %.3e', No, norm(d_Kt1,'fro'), n_Kt1));
    fs_expect(ok(d_Kt2,    n_Kt2), sprintf('node %d: Kt2 diff %.3e vs leg %.3e', No, norm(d_Kt2), n_Kt2));
    fs_expect(ok(d_Kn1(:), n_Kn1), sprintf('node %d: Kn1 diff %.3e vs leg %.3e', No, norm(d_Kn1,'fro'), n_Kn1));
    fs_expect(ok(d_Kn2,    n_Kn2), sprintf('node %d: Kn2 diff %.3e vs leg %.3e', No, norm(d_Kn2), n_Kn2));
end

fprintf('  [step 7] test done — calling fs_teardown for summary + exit code\n');
fs_teardown();
