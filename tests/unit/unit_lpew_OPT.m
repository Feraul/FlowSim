%UNIT_LPEW_OPT  Verify fs.lpew.OPT batched output equals legacy OPT_Interp_LPEW per-node slice.
%
%   Runs FlowSim's preprocessormod on M8.msh to get env, builds FS via
%   fs.mesh.build + fs.csr.buildCorners, then:
%     1. Calls fs.lpew.OPT once (batched, all nodes)
%     2. For a sample of nodes, calls legacy OPT_Interp_LPEW(No, env)
%     3. Asserts the batched slice equals the legacy output bit-exactly

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_lpew_OPT');

fsRoot = pwd;
setappdata(0, 'fs_test_root', fsRoot);
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

% Build a real env from M8.msh
fs_test_env('mesh', 'M8.msh');
env = preprocessormod(1);

% Adapt to FS
FS = fs.mesh.build(env);
FS = fs.csr.buildCorners(FS);

% Batched OPT
[P_all, T_all, O_all, Qc_all, Qn_all] = fs.lpew.OPT(FS);

% Structural checks
nNodes = FS.mesh.nNodes;
nCorners = FS.csr.nCorners;
nAllNeighbors = numel(FS.mesh.nsurn1);

fs_expect(size(P_all, 1) == nAllNeighbors, sprintf('P has %d rows (nAllNeighbors)', nAllNeighbors));
fs_expect(size(T_all, 1) == nAllNeighbors, 'T has nAllNeighbors rows');
fs_expect(size(O_all, 1) == nCorners,      sprintf('O has %d rows (nCorners)', nCorners));
fs_expect(size(Qc_all, 1) == nCorners,     'Qcorner has nCorners rows');
fs_expect(size(Qn_all, 1) == nAllNeighbors,'Qneighbor has nAllNeighbors rows');
fs_expect(size(P_all, 2) == 3 && size(O_all, 2) == 3, 'all outputs are 3-col');

% Per-node slice equality (sample interior + boundary nodes)
sampleNodes = [1, 5, 41, 45, 81];  % 4 corners of grid + one interior
for k = 1:numel(sampleNodes)
    No = sampleNodes(k);
    [O_leg, P_leg, T_leg, Qo_leg] = OPT_Interp_LPEW(No, env);

    % Slice from batched
    p_range = FS.mesh.nsurn2(No)+1 : FS.mesh.nsurn2(No+1);
    e_range = FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1);

    P_slice = P_all(p_range, :);
    T_slice = T_all(p_range, :);
    O_slice = O_all(e_range, :);
    Qc_slice = Qc_all(e_range(1), :);  % all rows are the same node's coord

    fs_expect(isequal(size(P_slice), size(P_leg)), sprintf('node %d: P size %dx%d matches legacy', No, size(P_slice,1), size(P_slice,2)));
    fs_expect(isequal(size(T_slice), size(T_leg)), sprintf('node %d: T size matches', No));
    fs_expect(isequal(size(O_slice), size(O_leg)), sprintf('node %d: O size matches', No));

    fs_expect(isequal(P_slice, P_leg), sprintf('node %d: P bit-identical', No));
    fs_expect(isequal(T_slice, T_leg), sprintf('node %d: T bit-identical', No));

    % O may have tiny float differences: precomputed centElem vs on-the-fly OPT sum
    % Compare via Frobenius relative diff (should be effectively zero — same formula)
    if isequal(size(O_slice), size(O_leg))
        rel = norm(O_slice - O_leg, 'fro') / max(1e-30, norm(O_leg, 'fro'));
        fs_expect(rel < 1e-12, sprintf('node %d: O Frobenius rel diff %.3e < 1e-12', No, rel));
    end

    fs_expect(isequal(Qc_slice, Qo_leg), sprintf('node %d: Qo bit-identical', No));
end

fs_teardown();
