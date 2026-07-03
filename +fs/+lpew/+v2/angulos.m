function [ve2, ve1, theta2, theta1] = angulos(FS, T_all, O_all, Q_corner)
%FS.LPEW.V2.ANGULOS  Batched LPEW2 corner angles (replaces angulos_Interp_LPEW2).
%
%   [ve2, ve1, theta2, theta1] = fs.lpew.v2.angulos(FS, T, O, Qcorner)
%
%   Inputs (from fs.lpew.OPT):
%     T          [nAllNeighbors x 3]  midpoints
%     O          [nCorners       x 3]  element centroids at each corner
%     Qcorner    [nCorners       x 3]  owning node's coord per corner
%
%   Preconditions (from fs.csr.buildCornerShifts):
%     FS.csr.tCurrent, FS.csr.tNext  [nCorners x 1]
%
%   Outputs (per-corner, one row per corner, flat):
%     ve1, ve2, theta1, theta2   each [nCorners x 1]
%
%   Per-node slice: legacy `angulos_Interp_LPEW2` returned 1×nec row vectors.
%   The equivalent slice is `ve1(FS.mesh.esurn2(No)+1 : FS.mesh.esurn2(No+1))`.

    tCurrent = FS.csr.tCurrent;
    tNext    = FS.csr.tNext;

    Tk  = T_all(tCurrent, :);
    Tk1 = T_all(tNext,    :);

    v0   = O_all - Q_corner;
    vth1 = Tk    - Q_corner;
    vth2 = Tk1   - Q_corner;
    v1   = Tk1   - Tk;

    n_v0   = sqrt(sum(v0.^2,   2));
    n_v1   = sqrt(sum(v1.^2,   2));
    n_vth1 = sqrt(sum(vth1.^2, 2));
    n_vth2 = sqrt(sum(vth2.^2, 2));

    d_vth1_v1 = -sum(vth1 .* v1, 2);   % dot(-vth1, v1)
    d_vth2_v1 =  sum(vth2 .* v1, 2);   % dot(-vth2, -v1) = dot(vth2, v1)
    d_v0_vth1 =  sum(v0   .* vth1, 2);
    d_v0_vth2 =  sum(v0   .* vth2, 2);

    ve1    = acos(d_vth1_v1 ./ (n_v1 .* n_vth1));
    ve2    = acos(d_vth2_v1 ./ (n_v1 .* n_vth2));
    theta1 = acos(d_v0_vth1 ./ (n_v0 .* n_vth1));
    theta2 = acos(d_v0_vth2 ./ (n_v0 .* n_vth2));
end
