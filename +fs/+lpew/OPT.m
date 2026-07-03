function [P, T, O, Qcorner, Qneighbor] = OPT(FS)
%FS.LPEW.OPT  Batched geometry gather for LPEW/interpolation.
%
%   [P, T, O, Qcorner, Qneighbor] = fs.lpew.OPT(FS)
%
%   Vectorized replacement for the per-node OPT_Interp_LPEW.m. Computes for
%   ALL nodes at once, in CSR-flat form.
%
%   Outputs:
%     P          [nAllNeighbors x 3]  neighbor-node coords
%     T          [nAllNeighbors x 3]  midpoints  (P + Qneighbor)/2
%     O          [nCorners       x 3]  element centroids at each corner
%     Qcorner    [nCorners       x 3]  owning node's coord per corner
%     Qneighbor  [nAllNeighbors  x 3]  owning node's coord per neighbor entry
%
%   Per-node slicing (for consumers that still want per-node views):
%     For node No:
%       P_No       = P(FS.mesh.nsurn2(No)+1 : FS.mesh.nsurn2(No+1), :)
%       T_No       = T(FS.mesh.nsurn2(No)+1 : FS.mesh.nsurn2(No+1), :)
%       O_No       = O(FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1), :)
%       Qo_No      = FS.mesh.coord(No, :)
%
%   Preconditions:
%     - FS.mesh.{coord, nsurn1, nsurn2, esurn1, esurn2} present
%     - FS.geom.centElem present (mesh-invariant, computed once by
%       preprocessormod; adapted via fs.mesh.build)
%     - FS.csr.cornerNode + cornerElem present (from fs.csr.buildCorners)
%
%   Complexity: O(nCorners + nAllNeighbors) indexed reads.
%   Eliminates the O(nNodes) loop present in legacy OPT_Interp_LPEW.

    coord      = FS.mesh.coord;
    nsurn1     = FS.mesh.nsurn1;
    nsurn2     = FS.mesh.nsurn2;
    centElem   = FS.geom.centElem;
    cornerNode = FS.csr.cornerNode;
    cornerElem = FS.csr.cornerElem;
    nNodes     = FS.mesh.nNodes;

    % ── Neighbors (P, T): one row per (node, neighbor) pair ──────────
    P = coord(nsurn1, :);
    nns = diff(nsurn2);
    nodeOfNeighbor = repelem((1:nNodes).', nns);
    Qneighbor = coord(nodeOfNeighbor, :);
    T = 0.5 * (P + Qneighbor);

    % ── Corners (O): one row per (node, element) pair ────────────────
    O = centElem(cornerElem, :);

    % Owning-node coord per corner (for downstream convenience)
    Qcorner = coord(cornerNode, :);
end
