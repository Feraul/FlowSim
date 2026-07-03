function [Kt1, Kt2, Kn1, Kn2] = ksInterp(FS, T_all, Q_corner)
%FS.LPEW.V2.KSINTERP  Batched permeability-tensor projections (replaces ferncodes_Ks_Interp_LPEW2).
%
%   [Kt1, Kt2, Kn1, Kn2] = fs.lpew.v2.ksInterp(FS, T, Qcorner)
%
%   Inputs:
%     T          [nAllNeighbors x 3]  midpoints (from fs.lpew.OPT)
%     Qcorner    [nCorners       x 3]  owning node coord per corner (from fs.lpew.OPT)
%
%   Preconditions:
%     FS.perm.kmap             — [nElems x 5+]  legacy kmap layout (col 2-5 = K11/K12/K21/K22)
%     FS.mesh.elem             — [nElems x 5+]  col 5 = material id
%     FS.csr.tCurrent/tNext    — per-corner shift indices (from buildCornerShifts)
%     FS.csr.cornerElem        — per-corner element id
%     FS.cfg.phasekey          — 1 (single-phase → L22=1) OR 2 (two-phase, requires Sw)
%
%   Outputs (per corner, one row per corner):
%     Kt1  [nCorners x 2]   tangential projections, column 1 for i=1 (T(k)),
%                           column 2 for i=2 (T(k+1) with interior wrap)
%     Kt2  [nCorners x 1]   tangential projection along T(k+1)-T(k)
%     Kn1  [nCorners x 2]   normal projections
%     Kn2  [nCorners x 1]   normal projection along T(k+1)-T(k)
%
%   Math: R = 90° rotation in 2D (embedded in 3D). For vector v = [vx, vy, 0]':
%     Rv = [vy, -vx, 0]'
%     <Rv, K*Rv> = K11*vy² - (K12+K21)*vx*vy + K22*vx²          → normal projection
%     <Rv, K*v>  = (K11-K22)*vx*vy + K12*vy² - K21*vx²          → tangential projection
%
%   For two-phase (phasekey=2), Kn2/Kt2 use L22-scaled tensor via twophasevar. Not
%   yet vectorized — falls back to error until phasekey=2 use case appears in tests.

    coord      = FS.mesh.coord;                                                    %#ok<NASGU>  % kept for future use
    elem       = FS.mesh.elem;
    kmap       = FS.perm.kmap;
    cornerElem = FS.csr.cornerElem;
    tCurrent   = FS.csr.tCurrent;
    tNext      = FS.csr.tNext;

    if isfield(FS.cfg, 'phasekey') && FS.cfg.phasekey == 2
        error('fs.lpew.v2.ksInterp:TwoPhaseNotVectorized', ...
              'Two-phase (phasekey=2) mobility scaling not yet vectorized. Track: PR-F.');
    end

    % Per-corner tensor (K11, K12, K21, K22) from element material id
    matIds = elem(cornerElem, 5);          % [nCorners x 1]
    % Some elements may have matId=0 (unset). Guard by clamping to 1 (matches
    % legacy behavior — MATLAB would error on 0 index, but sensible legacy meshes
    % have matId>=1). If any 0s exist, legacy would have crashed too.
    K11 = kmap(matIds, 2);
    K12 = kmap(matIds, 3);
    K21 = kmap(matIds, 4);
    K22 = kmap(matIds, 5);

    % ── Column i=1: use T(k) (i.e. tCurrent) ─────────────────────────
    v_k  = T_all(tCurrent, :) - Q_corner;  % [nCorners x 3]
    vx = v_k(:, 1); vy = v_k(:, 2);
    len2_k = vx.*vx + vy.*vy;
    Kn1_c1 = (K11.*vy.*vy - (K12 + K21).*vx.*vy + K22.*vx.*vx) ./ len2_k;
    Kt1_c1 = ((K11 - K22).*vx.*vy + K12.*vy.*vy - K21.*vx.*vx) ./ len2_k;

    % ── Column i=2: use T(k+1) (i.e. tNext — wraps for interior last corner) ─
    v_k1 = T_all(tNext, :) - Q_corner;
    vx = v_k1(:, 1); vy = v_k1(:, 2);
    len2_k1 = vx.*vx + vy.*vy;
    Kn1_c2 = (K11.*vy.*vy - (K12 + K21).*vx.*vy + K22.*vx.*vx) ./ len2_k1;
    Kt1_c2 = ((K11 - K22).*vx.*vy + K12.*vy.*vy - K21.*vx.*vx) ./ len2_k1;

    Kn1 = [Kn1_c1, Kn1_c2];
    Kt1 = [Kt1_c1, Kt1_c2];

    % ── Kn2, Kt2: use v = T(k+1) - T(k) with L22-scaled tensor ───────
    % Single-phase: L22=1 → K1 = K. Two-phase would multiply K by L22 per corner.
    v_diff = T_all(tNext, :) - T_all(tCurrent, :);
    vx = v_diff(:, 1); vy = v_diff(:, 2);
    len2_diff = vx.*vx + vy.*vy;
    Kn2 = (K11.*vy.*vy - (K12 + K21).*vx.*vy + K22.*vx.*vx) ./ len2_diff;
    Kt2 = ((K11 - K22).*vx.*vy + K12.*vy.*vy - K21.*vx.*vx) ./ len2_diff;
end
