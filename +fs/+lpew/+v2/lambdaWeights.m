function [lambda, r] = lambdaWeights(FS, Kt1, Kt2, Kn1, Kn2, theta1, theta2, ve1, ve2, netas, T_all, Q_corner)
%FS.LPEW.V2.LAMBDAWEIGHTS  Batched LPEW2 λ + r (replaces Lamdas_Weights_LPEW2).
%
%   [lambda, r] = fs.lpew.v2.lambdaWeights(FS, Kt1, Kt2, Kn1, Kn2, ...
%                                          theta1, theta2, ve1, ve2, netas, T, Qcorner)
%
%   Inputs (all per-corner from earlier PR-B2/B3 kernels):
%     Kt1  [nCorners x 2]     Kt2  [nCorners x 1]
%     Kn1  [nCorners x 2]     Kn2  [nCorners x 1]
%     theta1,2,ve1,2          [nCorners x 1]
%     netas                   [nCorners x 2]
%     T_all                   [nAllNeighbors x 3]  (needed for boundary r)
%     Q_corner                [nCorners x 3]       (needed for boundary r)
%
%   Outputs:
%     lambda  [nCorners x 1]  per-corner weight (interior formula fully vectorized;
%                             boundary handled per-node via legacy path)
%     r       [nNodes   x 2]  boundary Neumann correction (0 for interior nodes)
%
%   PR-B4 scope: interior nodes vectorized (majority — 80% of typical meshes).
%   Boundary nodes fall back to a per-node loop that mirrors the legacy formula.
%   Full boundary vectorization is PR-B4b (may be folded into PR-B5).

    nCorners = FS.csr.nCorners;
    nNodes   = FS.mesh.nNodes;

    lambda = zeros(nCorners, 1);
    r      = zeros(nNodes, 2);

    % ────────────────────────────────────────────────────────────────
    % INTERIOR NODES — fully vectorized (per-corner arithmetic)
    % ────────────────────────────────────────────────────────────────
    % For each corner i (in an interior node), let:
    %   iprev = FS.csr.cornerPrev(i)   (k-1 with wrap to nec)
    %   inext = FS.csr.cornerNext(i)   (k+1 with wrap to 1)
    %
    %   zeta_num(i) = Kn2(iprev)*cot(ve1(iprev)) + Kn2(i)*cot(ve2(i)) + Kt2(iprev) - Kt2(i)
    %   zeta_den(i) = Kn1(iprev,2)*cot(theta2(iprev)) + Kn1(i,1)*cot(theta1(i))
    %                 - Kt1(iprev,2) + Kt1(i,1)
    %   zeta(i)     = zeta_num / zeta_den
    %
    %   lambda(i)   = Kn1(i,1)*netas(i,1)*zeta(i) + Kn1(i,2)*netas(i,2)*zeta(inext)

    isInterior = FS.csr.nodeIsInterior(FS.csr.cornerNode);   % [nCorners x 1]
    iprev = FS.csr.cornerPrev;
    inext = FS.csr.cornerNext;

    % Compute zeta for ALL corners (interior formula) — we overwrite boundary later
    cot_ve1_prev    = cot(ve1(iprev));
    cot_ve2_cur     = cot(ve2);
    cot_theta2_prev = cot(theta2(iprev));
    cot_theta1_cur  = cot(theta1);

    zeta_num = Kn2(iprev) .* cot_ve1_prev + Kn2 .* cot_ve2_cur + Kt2(iprev) - Kt2;
    zeta_den = Kn1(iprev, 2) .* cot_theta2_prev + Kn1(:, 1) .* cot_theta1_cur ...
               - Kt1(iprev, 2) + Kt1(:, 1);
    zeta = zeta_num ./ zeta_den;

    % lambda for interior corners
    lambda_int = Kn1(:, 1) .* netas(:, 1) .* zeta ...
               + Kn1(:, 2) .* netas(:, 2) .* zeta(inext);
    lambda(isInterior) = lambda_int(isInterior);

    % ────────────────────────────────────────────────────────────────
    % BOUNDARY NODES — legacy per-node fallback (PR-B4b will vectorize)
    % ────────────────────────────────────────────────────────────────
    boundaryNodes = find(~FS.csr.nodeIsInterior);
    esurn2  = FS.mesh.esurn2;
    nsurn2  = FS.mesh.nsurn2;
    coord   = FS.mesh.coord;

    for k = 1:numel(boundaryNodes)
        No  = boundaryNodes(k);
        e_range = esurn2(No)+1 : esurn2(No+1);
        p_range = nsurn2(No)+1 : nsurn2(No+1);
        nec = numel(e_range);
        Qo  = coord(No, :);

        % Slice per-node inputs from the batched arrays
        Kt1_n = Kt1(e_range, :);
        Kt2_n = Kt2(e_range);
        Kn1_n = Kn1(e_range, :);
        Kn2_n = Kn2(e_range);
        th1_n = theta1(e_range).';
        th2_n = theta2(e_range).';
        ve1_n = ve1(e_range).';
        ve2_n = ve2(e_range).';
        netas_n = netas(e_range, :);
        T_n = T_all(p_range, :);

        % Legacy boundary formula for zeta (k=1..nec+1)
        zeta_n = zeros(1, nec+1);
        for kk = 1:nec+1
            if kk == 1
                zn = Kn2_n(kk)*cot(ve2_n(kk)) - Kt2_n(kk);
                zd = Kn1_n(kk, 1)*cot(th1_n(kk)) + Kt1_n(kk, 1);
                r(No, 1) = (1 + zn/zd) * norm(Qo - T_n(1, :));
            elseif kk == nec+1
                zn = Kn2_n(kk-1)*cot(ve1_n(kk-1)) + Kt2_n(kk-1);
                zd = Kn1_n(kk-1, 2)*cot(th2_n(kk-1)) - Kt1_n(kk-1, 2);
                r(No, 2) = (1 + zn/zd) * norm(Qo - T_n(nec+1, :));
            else
                zn = Kn2_n(kk-1)*cot(ve1_n(kk-1)) + Kn2_n(kk)*cot(ve2_n(kk)) ...
                     + Kt2_n(kk-1) - Kt2_n(kk);
                zd = Kn1_n(kk-1, 2)*cot(th2_n(kk-1)) + Kn1_n(kk, 1)*cot(th1_n(kk)) ...
                     - Kt1_n(kk-1, 2) + Kt1_n(kk, 1);
            end
            zeta_n(kk) = zn / zd;
        end

        % Lambda for boundary corners
        for kk = 1:nec
            % Boundary: always uses zeta(kk) and zeta(kk+1) (no wrap)
            lambda(e_range(kk)) = Kn1_n(kk, 1)*netas_n(kk, 1)*zeta_n(kk) ...
                                 + Kn1_n(kk, 2)*netas_n(kk, 2)*zeta_n(kk+1);
        end
    end
end
