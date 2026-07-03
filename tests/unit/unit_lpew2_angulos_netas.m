%UNIT_LPEW2_ANGULOS_NETAS  Verify fs.lpew.v2.{angulos,netas} vs legacy per-node.
%
%   Builds a real FS from M8.msh, calls the batched angulos + netas once,
%   then for a sample of nodes calls the legacy per-node versions and
%   asserts per-node slices match within 1e-12 relative Frobenius diff.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_lpew2_angulos_netas');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);

FS = fs.mesh.build(env);
FS = fs.csr.buildCorners(FS);
FS = fs.csr.buildCornerShifts(FS);

% One batched OPT call → gathers geometry for all nodes at once
[P_all, T_all, O_all, Qc_all, ~] = fs.lpew.OPT(FS);

% Batched angulos + netas
[ve2, ve1, theta2, theta1] = fs.lpew.v2.angulos(FS, T_all, O_all, Qc_all);
netas_all                  = fs.lpew.v2.netas(  FS, T_all, P_all, O_all, Qc_all);

fs_expect(numel(ve1) == FS.csr.nCorners, 've1 has nCorners rows');
fs_expect(numel(ve2) == FS.csr.nCorners, 've2 has nCorners rows');
fs_expect(numel(theta1) == FS.csr.nCorners, 'theta1 has nCorners rows');
fs_expect(numel(theta2) == FS.csr.nCorners, 'theta2 has nCorners rows');
fs_expect(isequal(size(netas_all), [FS.csr.nCorners, 2]), 'netas size correct');

% Verify per-node slice equality (sample: 4 corners + interior)
sampleNodes = [1, 5, 41, 45, 81];
tol = 1e-12;
for k = 1:numel(sampleNodes)
    No = sampleNodes(k);
    e_range = FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1);

    % Legacy per-node calls (need O, P, T, Qo from legacy OPT first)
    [O_leg, P_leg, T_leg, Qo_leg] = OPT_Interp_LPEW(No, env);
    [ve2_leg, ve1_leg, theta2_leg, theta1_leg] = ...
        angulos_Interp_LPEW2(O_leg, P_leg, T_leg, Qo_leg, No, env);
    netas_leg = netas_Interp_LPEW(O_leg, P_leg, T_leg, Qo_leg, No, env);

    % Batched slices — legacy returns row vectors (1×nec)
    ve1_slice    = ve1(e_range).';
    ve2_slice    = ve2(e_range).';
    theta1_slice = theta1(e_range).';
    theta2_slice = theta2(e_range).';
    netas_slice  = netas_all(e_range, :);

    fs_expect(isequal(size(ve1_slice), size(ve1_leg)), sprintf('node %d: ve1 size %dx%d', No, size(ve1_slice,1), size(ve1_slice,2)));

    rel_ve1    = norm(ve1_slice    - ve1_leg   ) / max(1e-30, norm(ve1_leg));
    rel_ve2    = norm(ve2_slice    - ve2_leg   ) / max(1e-30, norm(ve2_leg));
    rel_theta1 = norm(theta1_slice - theta1_leg) / max(1e-30, norm(theta1_leg));
    rel_theta2 = norm(theta2_slice - theta2_leg) / max(1e-30, norm(theta2_leg));
    rel_netas  = norm(netas_slice  - netas_leg,  'fro') / max(1e-30, norm(netas_leg, 'fro'));

    fs_expect(rel_ve1    < tol, sprintf('node %d: ve1 rel diff %.3e < %.0e',    No, rel_ve1,    tol));
    fs_expect(rel_ve2    < tol, sprintf('node %d: ve2 rel diff %.3e < %.0e',    No, rel_ve2,    tol));
    fs_expect(rel_theta1 < tol, sprintf('node %d: theta1 rel diff %.3e < %.0e', No, rel_theta1, tol));
    fs_expect(rel_theta2 < tol, sprintf('node %d: theta2 rel diff %.3e < %.0e', No, rel_theta2, tol));
    fs_expect(rel_netas  < tol, sprintf('node %d: netas rel diff %.3e < %.0e',  No, rel_netas,  tol));
end

fs_teardown();
