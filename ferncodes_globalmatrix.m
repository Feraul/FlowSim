%It is called by "ferncodes_solvepressure.m"

function [M,I,elembedge] = ferncodes_globalmatrix(env,preMPFAD,parmRichardEq)
% incializacao de parametros globais
auxcoord=env.geometry.coord;
auxelem=env.geometry.elem;
auxesurn2=env.geometry.esurn2;
auxesurn1=env.geometry.esurn1;
auxbedge=env.geometry.bedge;
auxinedge=env.geometry.inedge;
auxcentelem=env.geometry.centelem;
auxbcflag=env.config.bcflag;
auxnumcase=env.config.numcase;
auxnormals=env.geometry.normals;
auxmodflowcompared=env.config.modflowcase;
auxvicosity=env.config.visc;

% incializacao de parametros locais
Kt=preMPFAD.Kt;
Kn=preMPFAD.Kn;
Hesq=preMPFAD.Hesq;
Kde=preMPFAD.Kde;
Ded=preMPFAD.Ded;
%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(auxbedge,1);
inedgesize = size(auxinedge,1);

%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(auxelem,1),size(auxelem,1)); %Prealocação de M.
I = zeros(size(auxelem,1),1);
m=0;
jj=1;
elembedge=0;
% viscosidade ou mobilidade
visonface = 1;
%% ========================================================================
% -------------------------------------------------------------------------
% PRÉ-CÁLCULOS DE FLAGS E VISCOSIDADES
% -------------------------------------------------------------------------
isConc   = (200 < auxnumcase && auxnumcase < 300);
isSat    = (30  < auxnumcase && auxnumcase < 200);
is2phase = isConc && ismember(auxnumcase,[245 246 247 248 249 251]);

% boundary edges
v1b   = auxbedge(:,1);
v2b   = auxbedge(:,2);
elemL = auxbedge(:,3);
flagb = auxbedge(:,5);

v0_b  = auxcoord(v2b,:) - auxcoord(v1b,:);
v1_b  = auxcentelem(elemL,:) - auxcoord(v1b,:);
v2_b  = auxcentelem(elemL,:) - auxcoord(v2b,:);
nor_b = sqrt(sum((auxcoord(v1b,:) - auxcoord(v2b,:)).^2,2));

visonface_b = ones(bedgesize,1);
if isConc && is2phase
    visonface_b = auxviscosity(1:bedgesize,:);          % matriz (2 fases)
elseif isSat
    visonface_b = sum(auxviscosity(1:bedgesize,:),2);   % escalar por aresta
end

% internal edges
e1 = auxinedge(:,1);
e2 = auxinedge(:,2);
eL = auxinedge(:,3);
eR = auxinedge(:,4);

visonface_i = ones(inedgesize,1);
if isConc && is2phase
    visonface_i = auxviscosity(bedgesize + (1:inedgesize),:);
elseif isSat
    visonface_i = sum(auxviscosity(bedgesize + (1:inedgesize),:),2);
end

% -------------------------------------------------------------------------
% MONTAGEM VETORIZADA – BOUNDARY EDGES
% -------------------------------------------------------------------------
% Dirichlet vs Neumann
isDir = flagb < 200;
isNeu = ~isDir;

% ---------------- Dirichlet ----------------
% viscosidade efetiva (escalares)
if isConc && is2phase
    visD = sum(visonface_b(isDir,:),2);   % média simples (mantém lógica de uso escalar)
elseif isSat
    visD = visonface_b(isDir);
else
    visD = ones(nnz(isDir),1);
end

v0D   = v0_b(isDir,:);
v1D   = v1_b(isDir,:);
v2D   = v2_b(isDir,:);
norD  = nor_b(isDir);
lefD  = elemL(isDir);
v1n   = v1b(isDir);
v2n   = v2b(isDir);

c1D = preMPFAD.nflag(v1n,2);
c2D = preMPFAD.nflag(v2n,2);

KnD   = preMPFAD.Kn(isDir);
HesqD = preMPFAD.Hesq(isDir);
A_D   = -KnD ./ (HesqD .* sqrt(sum(v0D.^2,2)));

% contribuições em M(lef,lef)
rowsM_D = lefD;
colsM_D = lefD;
valsM_D = -visD .* A_D .* sum(v0D.^2,2);

% contribuições em I(lef)
dot_v2v0 = sum(v2D .* (-v0D),2);
dot_v1v0 = sum(v1D .*  v0D ,2);
valsI_D  = -visD .* A_D .* (dot_v2v0 .* c1D + dot_v1v0 .* c2D) + ...
    visD .* (c2D - c1D) .* Kt(isDir);

% ---------------- Neumann ----------------
lefN  = elemL(isNeu);
norN  = nor_b(isNeu);
flagN = flagb(isNeu);

valsI_N = zeros(nnz(isNeu),1);

if auxnumcase==341 || auxnumcase==341.1
    % ponto médio da aresta
    a1 = 0.5*(auxcoord(v1b(isNeu),:) + auxcoord(v2b(isNeu),:));
    if auxnumcase==341
        auxk = arrayfun(@(x,y) ferncodes_K(x,y), a1(:,1), a1(:,2));
    else
        auxk = arrayfun(@(x) ferncodes_K_1D(x), a1(:,1));
    end
    valsI_N = auxnormals(isNeu,2) .* auxk(:) .* nflagface(isNeu,2);
else
    % mapeia flagN em bcflag(:,1)
    [~,loc] = ismember(flagN, auxbcflag(:,1));
    valsI_N = norN .* auxbcflag(loc,2);
end

% -------------------------------------------------------------------------
% MONTAGEM VETORIZADA – INTERNAL EDGES
% -------------------------------------------------------------------------
% viscosidade efetiva interna
if isConc && is2phase
    visI = sum(visonface_i,2);
elseif isSat
    visI = visonface_i;
else
    visI = ones(inedgesize,1);
end

kI = visI .* preMPFAD.Kde;

% M contribuições (4 por aresta)
rowsM_I = [eL; eL; eR; eR];
colsM_I = [eL; eR; eR; eL];
valsM_I = [-kI; +kI; -kI; +kI];

% I contribuições por Dirichlet em vértices
maskD1 = preMPFAD.nflag(e1,1) < 200;
maskD2 = preMPFAD.nflag(e2,1) < 200;

valsI_L = zeros(size(I));
valsI_R = zeros(size(I));

% vértice 1
kD1 = kI(maskD1) .* preMPFAD.Ded(maskD1) .* preMPFAD.nflag(e1(maskD1),2);
idxL1 = eL(maskD1);
idxR1 = eR(maskD1);
valsI_L = accumarray(idxL1, -kD1, size(I));
valsI_R = valsI_R + accumarray(idxR1, +kD1, size(I));

% vértice 2
kD2 = kI(maskD2) .* preMPFAD.Ded(maskD2) .* preMPFAD.nflag(e2(maskD2),2);
idxL2 = eL(maskD2);
idxR2 = eR(maskD2);
valsI_L = valsI_L + accumarray(idxL2, +kD2, size(I));
valsI_R = valsI_R + accumarray(idxR2, -kD2, size(I));

% Neumann em vértices (201/202)
maskN1 = (preMPFAD.nflag(e1,1)==201 | preMPFAD.nflag(e1,1)==202);
maskN2 = (preMPFAD.nflag(e2,1)==201 | preMPFAD.nflag(e2,1)==202);

kN1 = kI(maskN1) .* preMPFAD.Ded(maskN1) .* preMPFAD.s(e1(maskN1));
idxL1N = eL(maskN1);
idxR1N = eR(maskN1);
valsI_L = valsI_L + accumarray(idxL1N, -kN1, size(I));
valsI_R = valsI_R + accumarray(idxR1N, +kN1, size(I));

kN2 = kI(maskN2) .* preMPFAD.Ded(maskN2) .* preMPFAD.s(e2(maskN2));
idxL2N = eL(maskN2);
idxR2N = eR(maskN2);
valsI_L = valsI_L + accumarray(idxL2N, +kN2, size(I));
valsI_R = valsI_R + accumarray(idxR2N, -kN2, size(I));

% -------------------------------------------------------------------------
% CONTRIBUIÇÕES DOS NÓS COM nflag > 200 (esurn1/esurn2)
% -------------------------------------------------------------------------
maskInt1 = preMPFAD.nflag(e1,1) > 200;
maskInt2 = preMPFAD.nflag(e2,1) > 200;

% vértice 1
edges1 = find(maskInt1);
rows_add = [];
cols_add = [];
vals_add = [];

for k = edges1.'
    n  = e1(k);
    idx = (auxesurn2(n)+1):auxesurn2(n+1);
    cols_loc = auxesurn1(idx);
    wloc     = preMPFAD.weight(idx).';
    rows_locL = repmat(eL(k),1,numel(cols_loc));
    rows_locR = repmat(eR(k),1,numel(cols_loc));
    valL = +visI(k)*preMPFAD.Kde(k)*preMPFAD.Ded(k)*wloc;
    valR = -visI(k)*preMPFAD.Kde(k)*preMPFAD.Ded(k)*wloc;
    rows_add = [rows_add, rows_locL, rows_locR];
    cols_add = [cols_add, cols_loc,  cols_loc ];
    vals_add = [vals_add, valL,      valR    ];
end

% vértice 2
edges2 = find(maskInt2);
for k = edges2.'
    n  = e2(k);
    idx = (auxesurn2(n)+1):auxesurn2(n+1);
    cols_loc = auxesurn1(idx);
    wloc     = preMPFAD.weight(idx).';
    rows_locL = repmat(eL(k),1,numel(cols_loc));
    rows_locR = repmat(eR(k),1,numel(cols_loc));
    valL = -visI(k)*preMPFAD.Kde(k)*preMPFAD.Ded(k)*wloc;
    valR = +visI(k)*preMPFAD.Kde(k)*preMPFAD.Ded(k)*wloc;
    rows_add = [rows_add, rows_locL, rows_locR];
    cols_add = [cols_add, cols_loc,  cols_loc ];
    vals_add = [vals_add, valL,      valR    ];
end

% -------------------------------------------------------------------------
% ASSEMBLAGEM FINAL EM M E I (ESPARSA)
% -------------------------------------------------------------------------
nelem = size(M,1);

% M: boundary + internal + esurn
rowsM_all = [rowsM_D; rowsM_I; rows_add(:)];
colsM_all = [colsM_D; colsM_I; cols_add(:)];
valsM_all = [valsM_D; valsM_I; vals_add(:)];

M = M + sparse(rowsM_all, colsM_all, valsM_all, nelem, nelem);

% I: boundary + internal
I = I + accumarray((1:nelem).',0, size(I)); % garante tamanho
I = I + accumarray(lefD, valsI_D, size(I));
I = I + accumarray(lefN, valsI_N, size(I));
I = I + valsI_L + valsI_R;

%==========================================================================
% calcula um problema transiente
% caso simulacao de aguas subterraneas, lei de Darcy
if (auxnumcase>330 || auxnumcase==330) && (auxnumcase<400)
    [M,I]=ferncodes_implicitandcranknicolson(M,I,env,dt);
end
% caso solo seco e fluido, Eq. Richards 
if 400<auxnumcase && auxnumcase<500
    [M,I]=soil_properties(M,I,parmRichardEq,preMPFAD,env);
end
%==========================================================================
% utilizase somente quando o teste vai ser comparado com resultados do modflow
if strcmp(auxmodflowcompared,'y')
    idx = elembedge(:,1);
    M(idx,:) = 0;
    M(sub2ind(size(M), idx, idx)) = 1;
    I(idx) = elembedge(:,2);
end
end
