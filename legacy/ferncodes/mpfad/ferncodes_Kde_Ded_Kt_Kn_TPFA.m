function [env] = ferncodes_Kde_Ded_Kt_Kn_TPFA(env, parms)

bedge    = env.geometry.bedge;
inedge   = env.geometry.inedge;
coord    = env.geometry.coord;
centelem = env.geometry.centelem;
elem     = env.geometry.elem;
nflagface= env.config.nflagface;

if env.config.numcase < 400 || isempty(parms)
    auxkmap = env.config.perm;
else
    auxkmap = parms.auxperm;
end

nb = size(bedge,1);
ni = size(inedge,1);
flowrateZ   = zeros(nb+ni,1);
flowresultZ = zeros(size(elem,1),1);

%% ── Boundary edges geometry ──────────────────────────────────────
v1     = coord(bedge(:,2),:) - coord(bedge(:,1),:);
Centro = centelem(bedge(:,3),:);
ve2    = coord(bedge(:,2),:) - Centro;

vx     = v1(:,1);
vy     = v1(:,2);
normv2 = vx.^2 + vy.^2;
nv     = sqrt(normv2);                              % substitui vecnorm

% cross 2D — so componente Z
Hesq   = abs(vx.*ve2(:,2) - vy.*ve2(:,1)) ./ nv;   % substitui cross+vecnorm

lef    = bedge(:,3);

%% ── Permeabilidade boundary — base ──────────────────────────────
matid = elem(bedge(:,3),5);
K11   = auxkmap(matid,2);
K12   = auxkmap(matid,3);
K21   = auxkmap(matid,4);
K22   = auxkmap(matid,5);

Kn = (K11.*vy.^2 - 2*K12.*vx.*vy + K22.*vx.^2) ./ normv2;

%% ── Flowrate boundary — delega ao benchmark ─────────────────────
if env.benchmark.temFlowrateBoundary()

    maskT      = bedge(:,5) < 200;
    h_contorno = nflagface(:,2);

    % ajusta K por modelo — substitui if numcase==437/438/else
    [K11, K12, K21, K22] = env.benchmark.ajustarKContorno(...
        env, parms, auxkmap, matid, h_contorno, maskT);

    Kn = (K11.*vy.^2 - 2*K12.*vx.*vy + K22.*vx.^2) ./ normv2;

    B1      = bedge(:,1);
    B2      = bedge(:,2);
    lef     = bedge(:,3);
    coordB1 = coord(B1,:);
    coordB2 = coord(B2,:);

    % norma da aresta — sem vecnorm
    dB      = coordB1 - coordB2;
    nor     = sqrt(dB(:,1).^2 + dB(:,2).^2);

    O       = centelem(lef,:);
    c1      = coord(B1,2);
    c2      = coord(B2,2);

    A1    = -Kn ./ (Hesq .* nv);
    term1 = sum((O-coordB2).*(coordB1-coordB2),2) .* c1;
    term2 = sum((O-coordB1).*(coordB2-coordB1),2) .* c2;

    % TPFA nao tem termo Kt
    flowrateZ(1:nb,1) = A1 .* (term1 + term2 - (nor.^2).*Centro(:,2));

    % Neumann
    mask201 = bedge(:,5) > 200;
    if any(mask201)
        flowrateZ(mask201) = flowrateZ(mask201);
    end

    % ajuste especifico do benchmark (ex: inversao caso 435)
    flowrateZ = env.benchmark.ajustarFlowrate(flowrateZ, bedge);

    flowresultZ = flowresultZ + accumarray(lef, flowrateZ(1:nb), size(flowresultZ));
end

%% ── Internal edges geometry ──────────────────────────────────────
C1   = centelem(inedge(:,3),:);
C2   = centelem(inedge(:,4),:);
lef  = inedge(:,3);
rel  = inedge(:,4);

vd1  = coord(inedge(:,2),:) - coord(inedge(:,1),:);
vd2  = C2 - coord(inedge(:,1),:);
ve2  = C1 - coord(inedge(:,1),:);

vx     = vd1(:,1);
vy     = vd1(:,2);
normv2 = vx.^2 + vy.^2;
nv     = sqrt(normv2);                              % substitui vecnorm

% cross 2D — so componente Z
H2 = abs(vx.*vd2(:,2) - vy.*vd2(:,1)) ./ nv;      % substitui cross+vecnorm
H1 = abs(vx.*ve2(:,2) - vy.*ve2(:,1)) ./ nv;

%% ── Permeabilidade internal ──────────────────────────────────────
matL = elem(inedge(:,3),5);
matR = elem(inedge(:,4),5);

K11L = auxkmap(matL,2);  K12L = auxkmap(matL,3);
K21L = auxkmap(matL,4);  K22L = auxkmap(matL,5);
K11R = auxkmap(matR,2);  K12R = auxkmap(matR,3);
K21R = auxkmap(matR,4);  K22R = auxkmap(matR,5);

%% ── Kn internal — TPFA nao usa Kt ───────────────────────────────
Kn1 = (K11L.*vy.^2 - 2*K12L.*vx.*vy + K22L.*vx.^2) ./ normv2;
Kn2 = (K11R.*vy.^2 - 2*K12R.*vx.*vy + K22R.*vx.^2) ./ normv2;

Kde = -nv .* (Kn1.*Kn2) ./ (Kn1.*H2 + Kn2.*H1);

%% ── Flowrate internal ────────────────────────────────────────────
if env.benchmark.temFlowrateBoundary()
    idx = nb + (1:ni);

    % TPFA: sem termo Ded (anisotropia nao capturada)
    flowrateZ(idx) = Kde .* (centelem(rel,2) - centelem(lef,2));

    flowresultZ = flowresultZ + ...
        accumarray(lef, flowrateZ(idx), size(flowresultZ)) - ...
        accumarray(rel, flowrateZ(idx), size(flowresultZ));
end

%% ── Empacota no env ──────────────────────────────────────────────
env.premethod.TPFA.Hesq        = Hesq;
env.premethod.TPFA.Kde         = Kde;
env.premethod.TPFA.Kn          = Kn;
env.premethod.TPFA.flowrateZ   = flowrateZ;
env.premethod.TPFA.flowresultZ = flowresultZ;
end