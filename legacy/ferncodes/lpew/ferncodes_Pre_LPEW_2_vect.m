function [env,weight,s] = ferncodes_Pre_LPEW_2_vect(env,parms)

if env.config.numcase > 400
    kmap = parms.auxperm;
else
    kmap = env.config.perm;
end

coord = env.geometry.coord;
elem  = env.geometry.elem;
ns1   = env.geometry.nsurn1;
ns2   = env.geometry.nsurn2;
es1   = env.geometry.esurn1;
es2   = env.geometry.esurn2;
N     = env.premethod.MPFAD.N;

nNodes = size(coord,1);
apw    = ones(nNodes+1,1);
r      = zeros(nNodes,2);
s      = zeros(nNodes,1);
weight = zeros(1,20*nNodes);

%% ── Pre-calcula centroide de cada elemento uma vez ───────────────
isQuad  = elem(:,4) ~= 0;
nvert   = 3 + isQuad;
centpre = zeros(size(elem,1),2);
for jj = 1:4
    mask = jj <= nvert;
    centpre(mask,1) = centpre(mask,1) + coord(elem(mask,jj),1)./nvert(mask);
    centpre(mask,2) = centpre(mask,2) + coord(elem(mask,jj),2)./nvert(mask);
end
sum_lambda_all = zeros(nNodes,1);
%% ── Loop sobre nos ───────────────────────────────────────────────
for y = 1:nNodes
    No  = y;
    Qox = coord(No,1);
    Qoy = coord(No,2);
    Qo  = [Qox Qoy 0];

    nns = ns2(No+1) - ns2(No);
    nec = es2(No+1) - es2(No);

    % ── P e T ─────────────────────────────────────────────────────
    P = zeros(nns,3);
    T = zeros(nns,3);
    ns_idx = ns2(No)+1 : ns2(No+1);
    P(:,1:2) = coord(ns1(ns_idx),1:2);
    T(:,1:2) = 0.5*(P(:,1:2) + Qo(1:2));

    % ── O — usa centpre pre-calculado ─────────────────────────────
    O = zeros(nec,3);
    es_idx = es2(No)+1 : es2(No+1);
    O(:,1:2) = centpre(es1(es_idx),:);

    %% ── Angulos ──────────────────────────────────────────────────
     ve2=zeros(1,es2(No+1)-es2(No));
    ve1=zeros(1,es2(No+1)-es2(No));
    theta2=zeros(1,es2(No+1)-es2(No));
    theta1=zeros(1,es2(No+1)-es2(No));


    for k=1:size(ve2,2),
        %Determinação dos vetores necessários à obtenção dos cossenos:
        v0=O(k,:)-Qo;
        if (k==size(ve2,2))&&(size(P,1)==size(O,1))
            vetorth2=T(1,:)-Qo;
            vetor1=T(1,:)-T(k,:);
        else
            vetor1=T(k+1,:)-T(k,:);
            vetorth2=T(k+1,:)-Qo;
        end
        vetorth1=T(k,:)-Qo;
        ve1(k)=acos(dot(-vetorth1,vetor1)/(norm(vetor1)*norm(vetorth1))); % revisar esses signos
        ve2(k)=acos(dot(-vetorth2,-vetor1)/(norm(vetor1)*norm(vetorth2)));
        theta2(k)=acos(dot(v0,vetorth2)/(norm(v0)*norm(vetorth2)));
        theta1(k)=acos(dot(v0,vetorth1)/(norm(v0)*norm(vetorth1)));
    end

    %% ── Netas — cross 2D (componente Z) sem cross() ─────────────
    netas = zeros(nec,2);

    if nns == nec   % no interno
        v1x = O(:,1) - Qox;   v1y = O(:,2) - Qoy;
        v2x = P(1:nec,1) - Qox; v2y = P(1:nec,2) - Qoy;

        % cross 2D: |v1 x v2| = |v1x*v2y - v1y*v2x|
        h1 = abs(v1x.*v2y - v1y.*v2x) ./ sqrt(v2x.^2 + v2y.^2);

        Tx = T(1:nec,1) - Qox;   Ty = T(1:nec,2) - Qoy;
        netas(:,1) = sqrt(Tx.^2 + Ty.^2) ./ h1;

        Psx = [P(2:nec,1); P(1,1)] - Qox;
        Psy = [P(2:nec,2); P(1,2)] - Qoy;
        h2  = abs(v1x.*Psy - v1y.*Psx) ./ sqrt(Psx.^2 + Psy.^2);

        Tsx = [T(2:nec,1); T(1,1)] - Qox;
        Tsy = [T(2:nec,2); T(1,2)] - Qoy;
        netas(:,2) = sqrt(Tsx.^2 + Tsy.^2) ./ h2;

    else   % no de contorno
        for k = 1:nec
            v1x = O(k,1)-Qox;  v1y = O(k,2)-Qoy;
            v2x = P(k,1)-Qox;  v2y = P(k,2)-Qoy;
            h1  = abs(v1x*v2y - v1y*v2x) / sqrt(v2x^2+v2y^2);
            Tkx = T(k,1)-Qox;  Tky = T(k,2)-Qoy;
            netas(k,1) = sqrt(Tkx^2+Tky^2) / h1;

            if k == nec
                v2x2 = P(end,1)-Qox; v2y2 = P(end,2)-Qoy;
                tvx  = T(end,1)-Qox; tvy  = T(end,2)-Qoy;
            else
                v2x2 = P(k+1,1)-Qox; v2y2 = P(k+1,2)-Qoy;
                tvx  = T(k+1,1)-Qox; tvy  = T(k+1,2)-Qoy;
            end
            h2 = abs(v1x*v2y2 - v1y*v2x2) / sqrt(v2x2^2+v2y2^2);
            netas(k,2) = sqrt(tvx^2+tvy^2) / h2;
        end
    end

    %% ── Tensores Kn1, Kt1, Kn2, Kt2 — R*v = [-vy, vx] ─────────
    % R = [0 1; -1 0] em 2D → R*[x;y] = [-y; x]
    Kn1 = zeros(nec,2);
    Kt1 = zeros(nec,2);
    Kn2 = zeros(nec,1);
    Kt2 = zeros(nec,1);

    % pre-indexa K para todos os elementos do no
    jj_all  = es1(es2(No)+1:es2(No+1));
    matids  = elem(jj_all,5);
    K11v = kmap(matids,2);
    K12v = kmap(matids,3);
    K21v = kmap(matids,4);
    K22v = kmap(matids,5);

    for k = 1:nec
        K11 = K11v(k); K12 = K12v(k); K21 = K21v(k); K22 = K22v(k);

        for i = 1:2
            if nns == nec && k == nec && i == 2
                Tv = T(1,1:2) - Qo(1:2);
            else
                Tv = T(k+i-1,1:2) - Qo(1:2);
            end
            % R*Tv = [-Tv(2), Tv(1)]
            Rx = -Tv(2);  Ry = Tv(1);
            n2 = Tv(1)^2 + Tv(2)^2;
            Kn1(k,i) = (Rx*(K11*Rx + K12*Ry) + Ry*(K21*Rx + K22*Ry)) / n2;
            Kt1(k,i) = (Rx*(K11*Tv(1)+K12*Tv(2)) + Ry*(K21*Tv(1)+K22*Tv(2))) / n2;
        end

        % Kn2, Kt2 — entre T(k) e T(k+1)
        if nns == nec && k == nec
            dT = T(1,1:2) - T(k,1:2);
        else
            dT = T(k+1,1:2) - T(k,1:2);
        end
        Rx  = -dT(2);  Ry = dT(1);
        n2  = dT(1)^2 + dT(2)^2;
        Kn2(k) = (Rx*(K11*Rx+K12*Ry) + Ry*(K21*Rx+K22*Ry)) / n2;
        Kt2(k) = (Rx*(K11*dT(1)+K12*dT(2)) + Ry*(K21*dT(1)+K22*dT(2))) / n2;
    end

    %% ── Lambda ───────────────────────────────────────────────────
    lambda = zeros(nec,1);
    zeta   = zeros(nec+1,1);

    if nns == nec   % no interno — shift vetorizado
        km1 = [nec, 1:nec-1];
        cve1 = cot(ve1);  cve2 = cot(ve2);
        cth1 = cot(theta1); cth2 = cot(theta2);

        zetan = Kn2(km1).*cve1(km1)' + Kn2.*cve2' + Kt2(km1) - Kt2;
        zetad = Kn1(km1,2).*cth2(km1)' + Kn1(:,1).*cth1' - Kt1(km1,2) + Kt1(:,1);
        zeta(1:nec) = zetan ./ zetad;

        kp1 = [2:nec, 1];
        lambda = zeta(1:nec).*Kn1(:,1).*netas(:,1) + ...
                 zeta(kp1)  .*Kn1(:,2).*netas(:,2);

    else   % no de contorno
        for k = 1:nec+1
            if k == 1
                zetan = Kn2(k)*cot(ve2(k)) - Kt2(k);
                zetad = Kn1(k,1)*cot(theta1(k)) + Kt1(k,1);
                r(No,1) = (1+zetan/zetad) * sqrt((Qox-T(1,1))^2+(Qoy-T(1,2))^2);
            elseif k == nec+1
                zetan = Kn2(k-1)*cot(ve1(k-1)) + Kt2(k-1);
                zetad = Kn1(k-1,2)*cot(theta2(k-1)) - Kt1(k-1,2);
                Tnp1x = T(min(nec+1,end),1); Tnp1y = T(min(nec+1,end),2);
                r(No,2) = (1+zetan/zetad) * sqrt((Qox-Tnp1x)^2+(Qoy-Tnp1y)^2);
            else
                zetan = Kn2(k-1)*cot(ve1(k-1)) + Kn2(k)*cot(ve2(k)) + Kt2(k-1) - Kt2(k);
                zetad = Kn1(k-1,2)*cot(theta2(k-1)) + Kn1(k,1)*cot(theta1(k)) ...
                      - Kt1(k-1,2) + Kt1(k,1);
            end
            zeta(k) = zetan/zetad;
        end
        for k = 1:nec
            if (k==nec)&&(size(P,1)==size(O,1))
                lambda(k)=zeta(k)*Kn1(k,1)*netas(k,1) + zeta(1)*Kn1(k,2)*netas(k,2);

            else
                lambda(k)=zeta(k)*Kn1(k,1)*netas(k,1) + zeta(k+1)*Kn1(k,2)*netas(k,2);

            end
        end
    end

    %% ── Pesos ────────────────────────────────────────────────────
    lambda = lambda(1:nec);
    sl     = sum(lambda);
    wloc   = lambda / sum(lambda);
    weight(1,apw(y):apw(y)+nec-1) = wloc;
    apw(y+1) = apw(y) + nec;
    % guarda sum_lambda 
    sum_lambda_all(y) = sl;
end
s = env.benchmark.calcularTermoNeumannVet(r, sum_lambda_all, N, env);
%% ── Saidas ───────────────────────────────────────────────────────
weight = weight(weight ~= 0);
env.premethod.MPFAD.weight = weight(1,1:apw(nNodes+1)-1);
env.premethod.MPFAD.s      = s;
end