function [M, I, elembedge] = build(env, parms)
%FS.ASSEMBLY.MPFAD.BUILD  Fully vectorized MPFA-D assembly.
%
%   [M, I, elembedge] = fs.assembly.mpfad.build(env, parms)
%
%   Drop-in replacement for ferncodes_globalmatrix_MPFAD. Legacy was already
%   ~90% vectorized; this version replaces the last 2 per-edge loops (LPEW
%   weight scatter over CSR corners) with fully-vectorized repelem gathers.

    coord      = env.geometry.coord;
    auxelem    = env.geometry.elem;
    auxesurn2  = env.geometry.esurn2;
    auxesurn1  = env.geometry.esurn1;
    bedge      = env.geometry.bedge;
    inedge     = env.geometry.inedge;
    auxcentelem= env.geometry.centelem;
    auxbcflag  = env.config.bcflag;
    auxnumcase = env.config.numcase;
    auxnormals = env.geometry.normals;
    auxmodflowcompared = env.config.modflowcase;
    nflag      = env.config.nflag;
    nflagface  = env.config.nflagface;

    pm         = env.premethod.MPFAD;
    Ded        = pm.Ded;
    flowresultZ= pm.flowresultZ;
    s          = pm.s;
    Kde        = pm.Kde;
    weight     = pm.weight;
    Kn         = pm.Kn;
    Hesq       = pm.Hesq;
    Kt         = pm.Kt;

    nelem      = size(auxelem, 1);
    bedgesize  = size(bedge, 1);
    inedgesize = size(inedge, 1);
    I          = zeros(nelem, 1);
    elembedge  = 0;

    isConc   = (200 < auxnumcase && auxnumcase < 300); %#ok<NASGU>
    isSat    = (30  < auxnumcase && auxnumcase < 200); %#ok<NASGU>

    % Boundary edges as arrays
    v1b   = bedge(:, 1);
    v2b   = bedge(:, 2);
    elemL = bedge(:, 3);
    flagb = bedge(:, 5);

    v0_b  = coord(v2b, :) - coord(v1b, :);
    v1_b  = auxcentelem(elemL, :) - coord(v1b, :);
    v2_b  = auxcentelem(elemL, :) - coord(v2b, :);
    nor_b = sqrt(sum((coord(v1b, :) - coord(v2b, :)).^2, 2));

    visonface_b = ones(bedgesize, 1); %#ok<NASGU>
    visonface_i = ones(inedgesize, 1); %#ok<NASGU>

    % ── Boundary — Dirichlet vs Neumann ─────────────────────────
    isDir = flagb < 200;
    isNeu = ~isDir;

    visD = ones(nnz(isDir), 1);
    v0D  = v0_b(isDir, :);
    v1D  = v1_b(isDir, :);
    v2D  = v2_b(isDir, :);
    lefD = elemL(isDir);
    v1n  = v1b(isDir);
    v2n  = v2b(isDir);

    c1D  = nflag(v1n, 2);
    c2D  = nflag(v2n, 2);
    KnD  = Kn(isDir);
    HesqD= Hesq(isDir);
    A_D  = -KnD ./ (HesqD .* sqrt(sum(v0D.^2, 2)));

    rowsM_D = lefD;
    colsM_D = lefD;
    valsM_D = -visD .* A_D .* sum(v0D.^2, 2);

    dot_v2v0 = sum(v2D .* (-v0D), 2);
    dot_v1v0 = sum(v1D .*  v0D , 2);
    valsI_D  = -visD .* A_D .* (dot_v2v0 .* c1D + dot_v1v0 .* c2D) + ...
               visD .* (c2D - c1D) .* Kt(isDir);

    lefN     = elemL(isNeu);
    valsI_N  = env.benchmark.calcularNeumannBoundary(isNeu, bedge, ...
        auxbcflag, nflagface, pm.flowrateZ, nor_b, auxnormals, env);

    % ── Internal edges ───────────────────────────────────────────
    e1 = inedge(:, 1);
    e2 = inedge(:, 2);
    eL = inedge(:, 3);
    eR = inedge(:, 4);

    visI = ones(inedgesize, 1);
    kI   = visI .* Kde;

    rowsM_I = [eL; eL; eR; eR];
    colsM_I = [eL; eR; eR; eL];
    valsM_I = [-kI;  kI; -kI;  kI];

    maskD1 = nflag(e1, 1) < 200;
    maskD2 = nflag(e2, 1) < 200;

    valsI_L = zeros(nelem, 1);
    valsI_R = zeros(nelem, 1);

    kD1 = kI(maskD1) .* Ded(maskD1) .* nflag(e1(maskD1), 2);
    valsI_L = valsI_L + accumarray(eL(maskD1), -kD1, [nelem 1]);
    valsI_R = valsI_R + accumarray(eR(maskD1),  kD1, [nelem 1]);

    kD2 = kI(maskD2) .* Ded(maskD2) .* nflag(e2(maskD2), 2);
    valsI_L = valsI_L + accumarray(eL(maskD2),  kD2, [nelem 1]);
    valsI_R = valsI_R + accumarray(eR(maskD2), -kD2, [nelem 1]);

    maskN1 = (nflag(e1, 1) == 201 | nflag(e1, 1) == 202);
    maskN2 = (nflag(e2, 1) == 201 | nflag(e2, 1) == 202);
    isBCvertex = (nflag(:, 1) == 201) | (nflag(:, 1) == 202);
    maskFaceHasBC  = isBCvertex(e1);
    idxFacesHasBC  = find(maskFaceHasBC);
    if ~isempty(idxFacesHasBC)
        kN1 = kI(idxFacesHasBC) .* Ded(idxFacesHasBC) .* s(e1(maskN1));
        idxL1N = inedge(idxFacesHasBC, 3);
        idxR1N = inedge(idxFacesHasBC, 4);
        valsI_L = valsI_L + accumarray(idxL1N, -kN1, [nelem 1]);
        valsI_R = valsI_R + accumarray(idxR1N,  kN1, [nelem 1]);
    end
    maskFaceHasBC1 = isBCvertex(e2);
    idxFacesHasBC1 = find(maskFaceHasBC1);
    if ~isempty(idxFacesHasBC1)
        kN2 = kI(idxFacesHasBC1) .* Ded(idxFacesHasBC1) .* s(e2(maskN2));
        idxL2N = inedge(idxFacesHasBC1, 3);
        idxR2N = inedge(idxFacesHasBC1, 4);
        valsI_L = valsI_L + accumarray(idxL2N,  kN2, [nelem 1]);
        valsI_R = valsI_R + accumarray(idxR2N, -kN2, [nelem 1]);
    end

    % ── FULLY VECTORIZED LPEW weight scatter ─────────────────────
    maskInt1 = nflag(e1, 1) > 200;
    maskInt2 = nflag(e2, 1) > 200;
    edges1   = find(maskInt1);
    edges2   = find(maskInt2);
    nec_per  = auxesurn2(2:end) - auxesurn2(1:end-1);

    [rows_add1, cols_add1, vals_add1] = weightScatter( ...
        edges1, e1, eL, eR, auxesurn1, auxesurn2, nec_per, weight, visI, Kde, Ded, +1);
    [rows_add2, cols_add2, vals_add2] = weightScatter( ...
        edges2, e2, eL, eR, auxesurn1, auxesurn2, nec_per, weight, visI, Kde, Ded, -1);

    % ── Final sparse assembly ────────────────────────────────────
    M = sparse( ...
        [rowsM_D; rowsM_I; rows_add1; rows_add2], ...
        [colsM_D; colsM_I; cols_add1; cols_add2], ...
        [valsM_D; valsM_I; vals_add1; vals_add2], ...
        nelem, nelem);

    I = I + accumarray(lefD, valsI_D, [nelem 1]);
    I = I + accumarray(lefN, valsI_N, [nelem 1]);
    I = I + valsI_L + valsI_R;

    [M, I] = env.benchmark.adicionarTermoTemporal(M, I, parms, flowresultZ, env);

    if strcmp(auxmodflowcompared, 'y')
        idx = elembedge(:, 1);
        M(idx, :) = 0;
        M(sub2ind(size(M), idx, idx)) = 1;
        I(idx) = elembedge(:, 2);
    end
end

function [rows_flat, cols_flat, vals_flat] = weightScatter( ...
    qedges, vArr, eL, eR, esurn1, esurn2, nec_per, weight, visI, Kde, Ded, signMul)
    if isempty(qedges)
        rows_flat = []; cols_flat = []; vals_flat = []; return;
    end
    vQ    = vArr(qedges);
    ncQ   = nec_per(vQ);
    startQ= esurn2(vQ) + 1;

    csum = cumsum([0; ncQ]);
    total = csum(end);
    if total == 0
        rows_flat = []; cols_flat = []; vals_flat = []; return;
    end
    posInGroup = (1:total).' - repelem(csum(1:end-1), ncQ) - 1;
    cornerFlat = repelem(startQ, ncQ) + posInGroup;

    cols_local = esurn1(cornerFlat);

    % weight is 1×nCorners row → transpose for column ops
    w_flat = weight(cornerFlat).';

    coef_per_edge = visI(qedges) .* Kde(qedges) .* Ded(qedges);
    coef_flat     = repelem(coef_per_edge, ncQ);

    rows_L = repelem(eL(qedges), ncQ);
    rows_R = repelem(eR(qedges), ncQ);
    % Legacy edges1 (signMul=+1): L = +coef*w, R = -coef*w
    % Legacy edges2 (signMul=-1): L = -coef*w, R = +coef*w
    vals_L = +signMul .* coef_flat .* w_flat;
    vals_R = -signMul .* coef_flat .* w_flat;

    rows_flat = [rows_L;  rows_R];
    cols_flat = [cols_local; cols_local];
    vals_flat = [vals_L;  vals_R];
end
