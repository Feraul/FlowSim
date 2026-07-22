%It is called by "ferncodes_solvepressure.m"

function [M,I,elembedge] = ferncodes_globalmatrix_MPFAD(env,parms)
% incializacao de parametros globais
coord=env.geometry.coord; elem=env.geometry.elem;esurn2=env.geometry.esurn2;
esurn1=env.geometry.esurn1;bedge=env.geometry.bedge;inedge=env.geometry.inedge;
centelem=env.geometry.centelem;bcflag=env.config.bcflag;
auxnumcase=env.config.numcase;normals=env.geometry.normals;
auxmodflowcompared=env.config.modflowcase;nflag=env.config.nflag;
Ded=env.premethod.MPFAD.Ded;flowrateZ=env.premethod.MPFAD.flowrateZ;
flowresultZ=env.premethod.MPFAD.flowresultZ;s=env.premethod.MPFAD.s;
Kde=env.premethod.MPFAD.Kde;weight=env.premethod.MPFAD.weight;
Kn=env.premethod.MPFAD.Kn;Hesq=env.premethod.MPFAD.Hesq;
nflagface = env.config.nflagface; Kt=env.premethod.MPFAD.Kt;
%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I = zeros(size(elem,1),1);
elembedge=0;
% viscosidade ou mobilidade
%% ========================================================================
% -------------------------------------------------------------------------
% PRÉ-CÁLCULOS DE FLAGS E VISCOSIDADES
% -------------------------------------------------------------------------
isConc   = (200 < auxnumcase && auxnumcase < 300);
isSat    = (30  < auxnumcase && auxnumcase < 200);
is2phase = isConc && ismember(auxnumcase,[245 246 247 248 249 251]);

% boundary edges
v1b   = bedge(:,1);
v2b   = bedge(:,2);
elemL = bedge(:,3);
flagb = bedge(:,5);

v0_b  = coord(v2b,:) - coord(v1b,:);
v1_b  = centelem(elemL,:) - coord(v1b,:);
v2_b  = centelem(elemL,:) - coord(v2b,:);
nor_b = sqrt(sum((coord(v1b,:) - coord(v2b,:)).^2,2));

visonface_b = ones(bedgesize,1);
if isConc && is2phase
    visonface_b = viscosity(1:bedgesize,:);          % matriz (2 fases)
elseif isSat
    visonface_b = sum(viscosity(1:bedgesize,:),2);   % escalar por aresta
end

% internal edges
e1 = inedge(:,1);
e2 = inedge(:,2);
eL = inedge(:,3);
eR = inedge(:,4);

visonface_i = ones(inedgesize,1);
if isConc && is2phase
    visonface_i = viscosity(bedgesize + (1:inedgesize),:);
elseif isSat
    visonface_i = sum(viscosity(bedgesize + (1:inedgesize),:),2);
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

c1D = nflag(v1n,2);
c2D = nflag(v2n,2);

KnD   = Kn(isDir);
HesqD = Hesq(isDir);
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

% ---------------- Neumann Boundary ----------------
lefN  = elemL(isNeu);

valsI_N = env.benchmark.calcularNeumannBoundary(isNeu, bedge, ...
    bcflag, nflagface, flowrateZ, nor_b, normals, env);

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

kI = visI .*Kde;

% M contribuições (4 por aresta)
rowsM_I = [eL; eL; eR; eR];
colsM_I = [eL; eR; eR; eL];
valsM_I = [-kI; +kI; -kI; +kI];

% I contribuições por Dirichlet em vértices
maskD1 = nflag(e1,1) < 200;
maskD2 = nflag(e2,1) < 200;

valsI_L = zeros(size(I));
valsI_R = zeros(size(I));

% vértice 1
kD1 = kI(maskD1) .* Ded(maskD1) .* nflag(e1(maskD1),2);
idxL1 = eL(maskD1);
idxR1 = eR(maskD1);
valsI_L = valsI_L+ accumarray(idxL1, -kD1, size(I));
valsI_R = valsI_R + accumarray(idxR1, +kD1, size(I));

% vértice 2
kD2 = kI(maskD2) .* Ded(maskD2) .* nflag(e2(maskD2),2);
idxL2 = eL(maskD2);
idxR2 = eR(maskD2);
valsI_L = valsI_L + accumarray(idxL2, +kD2, size(I));
valsI_R = valsI_R + accumarray(idxR2, -kD2, size(I));

% Neumann em vértices (201/202)
maskN1 = (nflag(e1,1)==201 | nflag(e1,1)==202);
maskN2 = (nflag(e2,1)==201 | nflag(e2,1)==202);
%--------------------------------------------------------------------------
% define quem é vértice de contorno (ajuste o critério aqui)
isBCvertex = (nflag(:,1) == 201) | (nflag(:,1) == 202);
% máscara de faces internas que têm pelo menos um vértice de contorno
maskFaceHasBC = isBCvertex(e1);

% índices das faces internas que tocam o contorno
idxFacesHasBC = find(maskFaceHasBC);
if ~isempty(idxFacesHasBC)
    %--------------------------------------------------------------------------
    kN1 = kI(idxFacesHasBC) .* Ded(idxFacesHasBC) .*s(e1(maskN1));
    idxL1N =inedge(idxFacesHasBC,3);
    idxR1N = inedge(idxFacesHasBC,4);
    valsI_L = valsI_L + accumarray(idxL1N, -kN1, size(I));
    valsI_R = valsI_R + accumarray(idxR1N, +kN1, size(I));
end
maskFaceHasBC1 =  isBCvertex(e2);

% índices das faces internas que tocam o contorno
idxFacesHasBC1 = find(maskFaceHasBC1);
if ~isempty(idxFacesHasBC1)
    kN2     = kI(idxFacesHasBC1) .* Ded(idxFacesHasBC1) .* s(e2(maskN2));
    idxL2N  = inedge(idxFacesHasBC1,3);
    idxR2N  = inedge(idxFacesHasBC1,4);
    valsI_L = valsI_L + accumarray(idxL2N, +kN2, size(I));
    valsI_R = valsI_R + accumarray(idxR2N, -kN2, size(I));
end

% -------------------------------------------------------------------------
% CONTRIBUIÇÕES DOS NÓS COM nflag > 200 (esurn1/esurn2), INCLUIDO NEUMANN
% -------------------------------------------------------------------------
maskInt1 = nflag(e1,1) > 200;
maskInt2 = nflag(e2,1) > 200;

% vértice 1
edges1 = find(maskInt1);
edges2= find(maskInt2);
rows_add = [];
cols_add = [];
vals_add = [];

% ── fora do loop — calcula tamanho total UMA vez ─────────────────
nec_per  = esurn2(2:end) - esurn2(1:end-1);
total    = 2*(sum(nec_per(e1(maskInt1))) + sum(nec_per(e2(maskInt2))));

% ── aloca UMA vez — sem realocação dentro do loop ─────────────────
rows_add = zeros(total,1);
cols_add = zeros(total,1);
vals_add = zeros(total,1);
ptr = 1;

for k = edges1.'
    n        = e1(k);
    idx      = esurn2(n)+1 : esurn2(n+1);
    nc       = numel(idx);
    cols_loc = esurn1(idx);
    coef     = visI(k)*Kde(k)*Ded(k);   % escalar — calcula uma vez

    % L — sem repmat, sem cópia extra
    ii = ptr:ptr+nc-1;
    rows_add(ii) = eL(k);               % escalar expandido automaticamente
    cols_add(ii) = cols_loc;
    vals_add(ii) = coef * weight(idx);  % sem wloc(:), sem transposta

    % R
    ii = ptr+nc : ptr+2*nc-1;
    rows_add(ii) = eR(k);
    cols_add(ii) = cols_loc;
    vals_add(ii) = -coef * weight(idx);

    ptr = ptr + 2*nc;
end

% vértice 2
% ── pre-aloca para edges2 ─────────────────────────────────────────
nec_per  = esurn2(2:end) - esurn2(1:end-1);
total2   = 2 * sum(nec_per(e2(maskInt2)));

rows_add2 = zeros(total2,1);
cols_add2 = zeros(total2,1);
vals_add2 = zeros(total2,1);
ptr2 = 1;

for k = edges2.'
    n        = e2(k);
    idx      = esurn2(n)+1 : esurn2(n+1);
    nc       = numel(idx);
    cols_loc = esurn1(idx);
    coef     = visI(k)*Kde(k)*Ded(k);   % escalar — calcula uma vez

    % L — sinal negativo para edges2
    ii = ptr2:ptr2+nc-1;
    rows_add2(ii) = eL(k);
    cols_add2(ii) = cols_loc;
    vals_add2(ii) = -coef * weight(idx);

    % R — sinal positivo para edges2
    ii = ptr2+nc : ptr2+2*nc-1;
    rows_add2(ii) = eR(k);
    cols_add2(ii) = cols_loc;
    vals_add2(ii) = +coef * weight(idx);

    ptr2 = ptr2 + 2*nc;
end

% -------------------------------------------------------------------------
% ASSEMBLAGEM FINAL EM M E I (ESPARSA)
% -------------------------------------------------------------------------
nelem = size(M,1);

% concatena os dois blocos pre-alocados
rows_all = [rows_add(1:ptr-1);  rows_add2(1:ptr2-1)];
cols_all = [cols_add(1:ptr-1);  cols_add2(1:ptr2-1)];
vals_all = [vals_add(1:ptr-1);  vals_add2(1:ptr2-1)];

% adiciona a M
M = M + sparse([rowsM_D; rowsM_I; rows_all], ...
               [colsM_D; colsM_I; cols_all], ...
               [valsM_D; valsM_I; vals_all], nelem, nelem);

% I: boundary + internal
I = I + accumarray((1:nelem).',0, size(I)); % garante tamanho
I = I + accumarray(lefD, valsI_D, size(I));
I = I + accumarray(lefN, valsI_N, size(I));
I = I + valsI_L + valsI_R;

%==========================================================================
% calcula um problema transiente
[M,I] = env.benchmark.adicionarTermoTemporal(M, I, parms, flowresultZ, env);
%==========================================================================
% utilizase somente quando o teste vai ser comparado com resultados do modflow
if strcmp(auxmodflowcompared,'y')
    idx = elembedge(:,1);
    M(idx,:) = 0;
    M(sub2ind(size(M), idx, idx)) = 1;
    I(idx) = elembedge(:,2);
end
end
