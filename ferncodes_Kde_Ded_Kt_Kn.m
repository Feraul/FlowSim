function [env] = ferncodes_Kde_Ded_Kt_Kn(env, parms)

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
B1  = bedge(:,1);  B2 = bedge(:,2);
v1x = coord(B2,1) - coord(B1,1);
v1y = coord(B2,2) - coord(B1,2);

Centro = centelem(bedge(:,3),:);
ve2x   = coord(B2,1) - Centro(:,1);
ve2y   = coord(B2,2) - Centro(:,2);

% cross 2D (componente Z) — substitui cross(v1,ve2,2) + vecnorm
normv2 = v1x.^2 + v1y.^2;
nv     = sqrt(normv2);                          % substitui vecnorm(v1,2,2)
Hesq   = abs(v1x.*ve2y - v1y.*ve2x) ./ nv;     % substitui vecnorm(ce,2,2)./nv

lef    = bedge(:,3);
vx     = v1x;
vy     = v1y;

%% ── Permeabilidade boundary ──────────────────────────────────────
matid = elem(bedge(:,3),5);
K11   = auxkmap(matid,2);
K12   = auxkmap(matid,3);
K21   = auxkmap(matid,4);
K22   = auxkmap(matid,5);

Kn = (K11.*vy.^2 - 2*K12.*vx.*vy + K22.*vx.^2) ./ normv2;
Kt = (vy.*(K11.*vx + K12.*vy) + (-vx).*(K21.*vx + K22.*vy)) ./ normv2;

%% ── Flowrate boundary ────────────────────────────────────────────
if env.benchmark.temFlowrateBoundary()

    maskT      = bedge(:,5) < 200;
    h_contorno = nflagface(:,2);

    [K11, K12, K21, K22] = env.benchmark.ajustarKContorno(...
        env, parms, auxkmap, matid, h_contorno, maskT);

    Kn = (K11.*vy.^2 - 2*K12.*vx.*vy + K22.*vx.^2) ./ normv2;
    Kt = (vy.*(K11.*vx + K12.*vy) + (-vx).*(K21.*vx + K22.*vy)) ./ normv2;

    coordB1 = coord(B1,:);
    coordB2 = coord(B2,:);
    dB      = coordB1 - coordB2;
    nor     = sqrt(dB(:,1).^2 + dB(:,2).^2);   % substitui sqrt(sum(...,2))

    O  = centelem(lef,:);
    c1 = coord(B1,2);
    c2 = coord(B2,2);

    A1    = -Kn ./ (Hesq .* nv);
    term1 = sum((O-coordB2).*(coordB1-coordB2),2) .* c1;
    term2 = sum((O-coordB1).*(coordB2-coordB1),2) .* c2;

    flowrateZ(1:nb) = A1.*(term1 + term2 - nor.^2.*Centro(:,2)) - (c2-c1).*Kt;

    % Neumann + ajuste especifico do benchmark
    mask201 = bedge(:,5) > 200;
    if any(mask201)
        flowrateZ(mask201) = flowrateZ(mask201);
    end
    flowrateZ = env.benchmark.ajustarFlowrate(flowrateZ, bedge);

    flowresultZ = flowresultZ + accumarray(lef, flowrateZ(1:nb), size(flowresultZ));
end

%% ── Internal edges geometry ──────────────────────────────────────
C1   = centelem(inedge(:,3),:);
C2   = centelem(inedge(:,4),:);
lef  = inedge(:,3);
rel  = inedge(:,4);
vcen = C2 - C1;

% vetores das arestas internas
e1x  = coord(inedge(:,2),1) - coord(inedge(:,1),1);
e1y  = coord(inedge(:,2),2) - coord(inedge(:,1),2);

vx     = e1x;
vy     = e1y;
normv2 = vx.^2 + vy.^2;
nv     = sqrt(normv2);                          % substitui vecnorm(vd1,2,2)

% vd2 = C2 - coord(inedge(:,1),:)
vd2x = C2(:,1) - coord(inedge(:,1),1);
vd2y = C2(:,2) - coord(inedge(:,1),2);

% ve2 = C1 - coord(inedge(:,1),:)
ve2x = C1(:,1) - coord(inedge(:,1),1);
ve2y = C1(:,2) - coord(inedge(:,1),2);

% cross 2D — substitui cross+vecnorm
H2 = abs(vx.*vd2y - vy.*vd2x) ./ nv;
H1 = abs(vx.*ve2y - vy.*ve2x) ./ nv;

no1 = coord(inedge(:,1),2);
no2 = coord(inedge(:,2),2);

%% ── Permeabilidade internal ──────────────────────────────────────
matL = elem(inedge(:,3),5);
matR = elem(inedge(:,4),5);

K11L = auxkmap(matL,2);  K12L = auxkmap(matL,3);
K21L = auxkmap(matL,4);  K22L = auxkmap(matL,5);
K11R = auxkmap(matR,2);  K12R = auxkmap(matR,3);
K21R = auxkmap(matR,4);  K22R = auxkmap(matR,5);

%% ── Kn Kt internal ───────────────────────────────────────────────
Kn1 = (K11L.*vy.^2 - 2*K12L.*vx.*vy + K22L.*vx.^2) ./ normv2;
Kt1 = (vy.*(K11L.*vx + K12L.*vy) + (-vx).*(K21L.*vx + K22L.*vy)) ./ normv2;
Kn2 = (K11R.*vy.^2 - 2*K12R.*vx.*vy + K22R.*vx.^2) ./ normv2;
Kt2 = (vy.*(K11R.*vx + K12R.*vy) + (-vx).*(K21R.*vx + K22R.*vy)) ./ normv2;

%% ── Transmissibilidades ──────────────────────────────────────────
Kde = -nv .* (Kn1.*Kn2) ./ (Kn1.*H2 + Kn2.*H1);

% dot(vd1,vcen,2) vetorizado sem dot()
dot_vd1_vcen = vx.*vcen(:,1) + vy.*vcen(:,2);
Ded = dot_vd1_vcen./normv2 - (1./nv).*((Kt2./Kn2).*H1 + (Kt1./Kn1).*H2);

%% ── Flowrate internal ────────────────────────────────────────────
if env.benchmark.temFlowrateBoundary()
    idx = nb + (1:ni);
    flowrateZ(idx) = Kde .* (centelem(rel,2) - centelem(lef,2) - Ded.*(no2-no1));
    flowresultZ    = flowresultZ + ...
        accumarray(lef, flowrateZ(idx), size(flowresultZ)) - ...
        accumarray(rel, flowrateZ(idx), size(flowresultZ));
end

%% ── Empacota ─────────────────────────────────────────────────────
env.premethod.MPFAD.Hesq        = Hesq;
env.premethod.MPFAD.Kde         = Kde;
env.premethod.MPFAD.Kn          = Kn;
env.premethod.MPFAD.Kt          = Kt;
env.premethod.MPFAD.Ded         = Ded;
env.premethod.MPFAD.flowrateZ   = flowrateZ;
env.premethod.MPFAD.flowresultZ = flowresultZ;
end