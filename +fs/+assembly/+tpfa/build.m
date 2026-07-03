function [M, I, elembedge] = build(env, parms)
%FS.ASSEMBLY.TPFA.BUILD  Vectorized TPFA assembly (drop-in for ferncodes_globalmatrix_TPFA).
    coord      = env.geometry.coord;
    auxelem    = env.geometry.elem;
    bedge      = env.geometry.bedge;
    inedge     = env.geometry.inedge;
    auxcentelem= env.geometry.centelem;
    auxbcflag  = env.config.bcflag;
    auxnumcase = env.config.numcase;
    auxnormals = env.geometry.normals;
    auxmodflowcompared = env.config.modflowcase;
    nflag      = env.config.nflag;
    nflagface  = env.config.nflagface;
    flowrateZ  = env.premethod.TPFA.flowrateZ;
    flowresultZ= env.premethod.TPFA.flowresultZ;

    nelem      = size(auxelem, 1);
    inedgesize = size(inedge, 1);
    I          = zeros(nelem, 1);
    elembedge  = 0;
    dt         = 1;

    v1b   = bedge(:, 1);
    v2b   = bedge(:, 2);
    elemL = bedge(:, 3);
    flagb = bedge(:, 5);

    v0_b  = coord(v2b, :) - coord(v1b, :);
    v1_b  = auxcentelem(elemL, :) - coord(v1b, :);
    v2_b  = auxcentelem(elemL, :) - coord(v2b, :);
    nor_b = sqrt(sum((coord(v1b, :) - coord(v2b, :)).^2, 2));

    eL = inedge(:, 3);
    eR = inedge(:, 4);

    isDir = flagb < 200;
    isNeu = ~isDir;

    visD  = ones(nnz(isDir), 1);
    v0D   = v0_b(isDir, :);
    v1D   = v1_b(isDir, :);
    v2D   = v2_b(isDir, :);
    lefD  = elemL(isDir);
    v1n   = v1b(isDir);
    v2n   = v2b(isDir);
    c1D   = nflag(v1n, 2);
    c2D   = nflag(v2n, 2);
    KnD   = env.premethod.TPFA.Kn(isDir);
    HesqD = env.premethod.TPFA.Hesq(isDir);
    A_D   = -KnD ./ (HesqD .* sqrt(sum(v0D.^2, 2)));

    rowsM_D = lefD;
    colsM_D = lefD;
    valsM_D = -visD .* A_D .* sum(v0D.^2, 2);
    dot_v2v0 = sum(v2D .* (-v0D), 2);
    dot_v1v0 = sum(v1D .*  v0D , 2);
    valsI_D  = -visD .* A_D .* (dot_v2v0 .* c1D + dot_v1v0 .* c2D);

    lefN  = elemL(isNeu);
    norN  = nor_b(isNeu);
    flagN = flagb(isNeu);
    valsI_N = zeros(nnz(isNeu), 1);
    if auxnumcase == 341 || auxnumcase == 341.1
        a1 = 0.5 * (coord(v1b(isNeu), :) + coord(v2b(isNeu), :));
        if auxnumcase == 341
            auxk = arrayfun(@(x, y) ferncodes_K(x, y), a1(:, 1), a1(:, 2));
        else
            auxk = arrayfun(@(x) ferncodes_K_1D(x), a1(:, 1));
        end
        valsI_N = auxnormals(isNeu, 2) .* auxk(:) .* nflagface(isNeu, 2);
    else
        [~, loc] = ismember(flagN, auxbcflag(:, 1));
        mask222 = bedge(:, 5) > 200;
        valsI_N = norN .* auxbcflag(loc, 2) + flowrateZ(find(mask222 == 1), 1);
    end

    visI = ones(inedgesize, 1);
    kI   = visI .* env.premethod.TPFA.Kde;
    rowsM_I = [eL; eL; eR; eR];
    colsM_I = [eL; eR; eR; eL];
    valsM_I = [-kI;  kI; -kI;  kI];

    M = sparse( ...
        [rowsM_D; rowsM_I], ...
        [colsM_D; colsM_I], ...
        [valsM_D; valsM_I], nelem, nelem);
    I = I + accumarray(lefD, valsI_D, [nelem 1]);
    I = I + accumarray(lefN, valsI_N, [nelem 1]);

    if (auxnumcase > 330 || auxnumcase == 330) && (auxnumcase < 400)
        [M, I] = ferncodes_implicitandcranknicolson(M, I, env, dt);
    end
    if 400 < auxnumcase && auxnumcase < 500
        [M, I] = soil_properties(M, I, parms, flowresultZ, env);
    end

    if strcmp(auxmodflowcompared, 'y')
        idx = elembedge(:, 1);
        M(idx, :) = 0;
        M(sub2ind(size(M), idx, idx)) = 1;
        I(idx) = elembedge(:, 2);
    end
end
