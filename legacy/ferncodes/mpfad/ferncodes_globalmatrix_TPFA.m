
%--------------------------------------------------------------------------

function [M,I,elembedge] =ferncodes_globalmatrix_TPFA(env,parms)
% incializacao de parametros globais
coord=env.geometry.coord;
auxelem=env.geometry.elem;

bedge=env.geometry.bedge;
inedge=env.geometry.inedge;
centelem=env.geometry.centelem;
bcflag=env.config.bcflag;
numcase=env.config.numcase;
normals=env.geometry.normals;
modflowcompared=env.config.modflowcase;
nflag=env.config.nflag;
nflagface=env.config.nflagface;
flowrateZ=env.premethod.TPFA.flowrateZ;
flowresultZ=env.premethod.TPFA.flowresultZ;


%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

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
isConc   = (200 < numcase && numcase < 300);
isSat    = (30  < numcase && numcase < 200);
is2phase = isConc && ismember(numcase,[245 246 247 248 249 251]);

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
    visonface_b = auxviscosity(1:bedgesize,:);          % matriz (2 fases)
elseif isSat
    visonface_b = sum(auxviscosity(1:bedgesize,:),2);   % escalar por aresta
end

% internal edges
e1 = inedge(:,1);
e2 = inedge(:,2);
eL = inedge(:,3);
eR = inedge(:,4);

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

c1D = nflag(v1n,2);
c2D = nflag(v2n,2);

KnD   = env.premethod.TPFA.Kn(isDir);
HesqD = env.premethod.TPFA.Hesq(isDir);
A_D   = -KnD ./ (HesqD .* sqrt(sum(v0D.^2,2)));

% contribuições em M(lef,lef)
rowsM_D = lefD;
colsM_D = lefD;
valsM_D = -visD .* A_D .* sum(v0D.^2,2);

% contribuições em I(lef)
dot_v2v0 = sum(v2D .* (-v0D),2);
dot_v1v0 = sum(v1D .*  v0D ,2);
valsI_D  = -visD .* A_D .* (dot_v2v0 .* c1D + dot_v1v0 .* c2D);

% ---------------- Neumann ----------------
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

kI = visI .* env.premethod.TPFA.Kde;

% M contribuições (4 por aresta)
rowsM_I = [eL; eL; eR; eR];
colsM_I = [eL; eR; eR; eL];
valsM_I = [-kI; +kI; -kI; +kI];

% -------------------------------------------------------------------------
% ASSEMBLAGEM FINAL EM M E I (ESPARSA)
% -------------------------------------------------------------------------
nelem = size(M,1);

% M: boundary + internal + esurn
rowsM_all = [rowsM_D; rowsM_I];
colsM_all = [colsM_D; colsM_I];
valsM_all = [valsM_D; valsM_I];

M = M + sparse(rowsM_all, colsM_all, valsM_all, nelem, nelem);

% I: boundary + internal
I = I + accumarray((1:nelem).',0, size(I)); % garante tamanho
I = I + accumarray(lefD, valsI_D, size(I));
I = I + accumarray(lefN, valsI_N, size(I));

%==========================================================================
%==========================================================================
% calcula um problema transiente
[M,I] = env.benchmark.adicionarTermoTemporal(M, I, parms, flowresultZ, env);
%==========================================================================
% utilizase somente quando o teste vai ser comparado com resultados do modflow
if strcmp(modflowcompared,'y')
    idx = elembedge(:,1);
    M(idx,:) = 0;
    M(sub2ind(size(M), idx, idx)) = 1;
    I(idx) = elembedge(:,2);
end
end
