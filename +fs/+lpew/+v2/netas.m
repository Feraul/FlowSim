function out = netas(FS, T_all, P_all, O_all, Q_corner)
%FS.LPEW.V2.NETAS  Batched LPEW netas (replaces netas_Interp_LPEW).
%
%   netas = fs.lpew.v2.netas(FS, T, P, O, Qcorner)   → [nCorners x 2]
%
%   Inputs (from fs.lpew.OPT):
%     T, P        [nAllNeighbors x 3]
%     O           [nCorners       x 3]
%     Qcorner     [nCorners       x 3]
%
%   Preconditions (from fs.csr.buildCornerShifts):
%     FS.csr.tCurrent, FS.csr.tNext

    tCurrent = FS.csr.tCurrent;
    tNext    = FS.csr.tNext;

    Tk  = T_all(tCurrent, :);
    Tk1 = T_all(tNext,    :);
    Pk  = P_all(tCurrent, :);
    Pk1 = P_all(tNext,    :);

    v1      = O_all - Q_corner;
    v2_col1 = Pk    - Q_corner;
    v2_col2 = Pk1   - Q_corner;

    ce1 = cross(v1, v2_col1, 2);
    ce2 = cross(v1, v2_col2, 2);

    n_ce1     = sqrt(sum(ce1.^2,     2));
    n_ce2     = sqrt(sum(ce2.^2,     2));
    n_v2_col1 = sqrt(sum(v2_col1.^2, 2));
    n_v2_col2 = sqrt(sum(v2_col2.^2, 2));

    h1 = n_ce1 ./ n_v2_col1;
    h2 = n_ce2 ./ n_v2_col2;

    Tk_dist  = sqrt(sum((Tk  - Q_corner).^2, 2));
    Tk1_dist = sqrt(sum((Tk1 - Q_corner).^2, 2));

    out = [Tk_dist ./ h1, Tk1_dist ./ h2];
end
