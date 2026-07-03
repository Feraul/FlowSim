%UNIT_CSR_CORNERS  Verify fs.csr.buildCorners produces correct CSR-flat layout.
addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_csr_corners');

% ── Case: 2×2 quad mesh (fabricated). 9 nodes, 4 elements ─────────
% Node numbering:   7 8 9
%                   4 5 6
%                   1 2 3
% Elements (CCW): 1={1,2,5,4}  2={2,3,6,5}  3={4,5,8,7}  4={5,6,9,8}
%
% esurn (elements-around-node), CSR form:
% node 1: {1}          → 1 corner
% node 2: {1,2}        → 2 corners
% node 3: {2}          → 1 corner
% node 4: {1,3}        → 2 corners
% node 5: {1,2,3,4}    → 4 corners
% node 6: {2,4}        → 2 corners
% node 7: {3}          → 1 corner
% node 8: {3,4}        → 2 corners
% node 9: {4}          → 1 corner
%
% Total corners = 1+2+1+2+4+2+1+2+1 = 16
esurn1 = [1; 1;2; 2; 1;3; 1;2;3;4; 2;4; 3; 3;4; 4];
esurn2 = [0; 1; 3; 4; 6; 10; 12; 13; 15; 16];  % length nNodes+1 = 10

FS.mesh = struct('nNodes', 9, 'nElems', 4, 'esurn1', esurn1, 'esurn2', esurn2);
FS.csr  = struct();

FS = fs.csr.buildCorners(FS);

fs_expect(FS.csr.nCorners == 16,           'nCorners = 16');
fs_expect(numel(FS.csr.nodePtr) == 10,     'nodePtr length = nNodes+1');
fs_expect(numel(FS.csr.cornerNode) == 16,  'cornerNode length = nCorners');
fs_expect(numel(FS.csr.cornerElem) == 16,  'cornerElem length = nCorners');
fs_expect(isequal(FS.csr.cornerElem, esurn1), 'cornerElem = esurn1');
fs_expect(FS.csr.maxNec == 4,              'maxNec = 4 (interior node has 4 corners)');
fs_expect(sum(FS.csr.nodeNec) == 16,       'sum(nodeNec) = nCorners');

% cornerNode invariant: element k belongs to node cornerNode(k),
% and the sequence of cornerNode should be nondecreasing (CSR order).
fs_expect(all(diff(FS.csr.cornerNode) >= 0), 'cornerNode is nondecreasing (CSR order)');

% Each node i appears exactly nodeNec(i) times in cornerNode
counts = accumarray(FS.csr.cornerNode, 1, [FS.mesh.nNodes 1]);
fs_expect(isequal(counts, FS.csr.nodeNec), 'accumarray(cornerNode) reproduces nodeNec');

% Interior node 5 has 4 corners, all in elements {1,2,3,4}
mask5 = FS.csr.cornerNode == 5;
fs_expect(isequal(sort(FS.csr.cornerElem(mask5)), [1;2;3;4]), 'node 5 corners cover elements {1,2,3,4}');

fs_teardown();
