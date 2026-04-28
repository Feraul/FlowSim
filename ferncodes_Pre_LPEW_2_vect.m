function [preMPFAD,weight,s] = ferncodes_Pre_LPEW_2_vect(preMPFAD,parmRichardEq,env)
auxnumcase = env.config.numcase;

if auxnumcase > 400
    kmap = parmRichardEq.auxperm;
else
    kmap = env.config.perm;
end

coord = env.geometry.coord;
elem  = env.geometry.elem;

ns1 = env.geometry.nsurn1;
ns2 = env.geometry.nsurn2;

es1 = env.geometry.esurn1;
es2 = env.geometry.esurn2;
N=preMPFAD.N;

nNodes = size(coord,1);

apw = ones(nNodes+1,1);
r   = zeros(nNodes,2);
s   = 0;

weight = zeros(1,20*nNodes); % tamanho seguro
R = [0 1 0; -1 0 0; 0 0 0];
% vetor de pontos na vizinhança do nó "ni".

%% -------------------------------------------------
% Loop sobre nós
%% -------------------------------------------------
for y = 1:nNodes
    No = y;
    Qo = coord(No,:);
    P=zeros(ns2(No+1)-ns2(No),3);
    T=zeros(ns2(No+1)-ns2(No),3);
    O=zeros(es2(No+1)-es2(No),3);
    nec=size(O,1);
    for i=1:size(P,1),
        P(i,:)=coord(ns1(ns2(No)+i),:);
        T(i,:)=(P(i,:)+Qo)/2;
        % aloca o flag da da aresta ou face pertence ao contorno
    end

    %Construção do vetor O, dos centróides (pontos de colocação) dos elementos%
    %que concorrem no nó ni.                                                  %

    for i=1:size(O,1),
        %Verifica se o elemento é um quadrilátero ou um triângulo.
        if elem(es1(es2(No)+i),4)==0 % lenbrando que o quarta columna
            b=3;
        else
            b=4;  % da matriz de elementos é para quadrilateros
        end
        %Carrega adequadamente o vetor O (braicentro de cada elemento)
        for j=1:b
            O(i,1)=O(i,1)+(coord(elem(es1(es2(No)+i),j),1)/b);
            O(i,2)=O(i,2)+(coord(elem(es1(es2(No)+i),j),2)/b);
            O(i,3)=O(i,3)+(coord(elem(es1(es2(No)+i),j),3)/b);
        end
    end
    %% -------------------------------------------------
    % Vetores para ângulos
    %% -------------------------------------------------
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

    %% -------------------------------------------------
    % Vetores netas
    %% -------------------------------------------------
    netas = zeros(nec,2);

    if size(P,1) == size(O,1)
        % Nó interno (vectorizado)
        v1 = O - Qo;
        v2 = P - Qo;
        ce = cross(v1,v2,2);
        h1 = vecnorm(ce,2,2)./vecnorm(v2,2,2);
        netas(:,1) = vecnorm(T-Qo,2,2)./h1;

        Pshift = [P(2:end,:); P(1,:)];
        Tshift = [T(2:end,:); T(1,:)];
        v2 = Pshift - Qo;
        ce = cross(v1,v2,2);
        h2 = vecnorm(ce,2,2)./vecnorm(v2,2,2);
        netas(:,2) = vecnorm(Tshift-Qo,2,2)./h2;
    else
        % Nó de contorno (loop)
        for k = 1:nec
            v1c = O(k,:) - Qo;
            v2c = P(k,:) - Qo;
            ce = cross(v1c,v2c);
            h1 = norm(ce)/norm(v2c);
            netas(k,1) = norm(T(k,:) - Qo)/h1;

            if k==nec
                v2c = P(end,:) - Qo;
                tvec = T(end,:) - Qo;
            else
                v2c = P(k+1,:) - Qo;
                tvec = T(k+1,:) - Qo;
            end
            ce = cross(v1c,v2c);
            h2 = norm(ce)/norm(v2c);
            netas(k,2) = norm(tvec)/h2;
        end
    end

    %% -------------------------------------------------
    % Tensores Kn1, Kt1, Kn2, Kt2
    %% -------------------------------------------------

    Kt1=zeros(nec,2); %As colunas representam i=1 e i=2.
    Kt2=zeros(nec,1);
    Kn1=zeros(nec,2);
    Kn2=zeros(nec,1);
    K=zeros(3);
    R=[0 1 0; -1 0 0; 0 0 0];
    %Construção do tensor permeabilidade.%

    %Cálculo das primeiras constantes, para todas as células que concorrem num
    %vertice "No".
    for k=1:nec
        % elemento j
        j=es1(es2(No)+k);
        % permeabilidade
        K(1,1)=kmap(elem(j,5),2);
        K(1,2)=kmap(elem(j,5),3);
        K(2,1)=kmap(elem(j,5),4);
        K(2,2)=kmap(elem(j,5),5);
        for i=1:2
            if (size(T,1)==size(O,1))&&(k==nec)&&(i==2)
                Kn1(k,i)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/norm(T(1,:)-Qo)^2;
                Kt1(k,i)=((R*(T(1,:)-Qo)')'*K*(T(1,:)-Qo)')/norm(T(1,:)-Qo)^2;
            else
                Kn1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(R*(T(k+i-1,:)-Qo)'))/norm(T(k+i-1,:)-Qo)^2;
                Kt1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(T(k+i-1,:)-Qo)')/norm(T(k+i-1,:)-Qo)^2;
            end
        end
        %------------------------- Tensores ----------------------------------%
        if (size(T,1)==size(O,1))&&(k==nec)
            %------------ Calculo dos K's internos no elemento ---------------%
            Kn2(k)=((R*(T(1,:)-T(k,:))')'*K*(R*(T(1,:)-T(k,:))'))/norm(T(1,:)-T(k,:))^2;
            Kt2(k)=((R*(T(1,:)-T(k,:))')'*K*(T(1,:)-T(k,:))')/norm(T(1,:)-T(k,:))^2;

        else
            Kn2(k)=(R*(T(k+1,:)-T(k,:))')'*K*(R*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
            Kt2(k)=((R*(T(k+1,:)-T(k,:))')'*K*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;

        end
    end

    %% -------------------------------------------------
    % Cálculo de lambda
    %% -------------------------------------------------
    %Determina os lambdas.
    nec=size(O,1);
    lambda=zeros(nec,1);
    % see eq 26 of the article: "A linearity-preserving cell-centered scheme for the heterogeneous
    %and anisotropic diffusion equations on general meshes" DOI: 10.1002/fld
    if size(P,1)==size(O,1) %Se for um nó interno.
        for k=1:nec,
            if (k==1)&&(size(P,1)==size(O,1))
                zetan=Kn2(nec)*cot(ve1(nec))+Kn2(k)*cot(ve2(k))+Kt2(nec)-Kt2(k);
                zetad=Kn1(nec,2)*cot(theta2(nec))+Kn1(k,1)*cot(theta1(k)) ...
                    -Kt1(nec,2)+Kt1(k,1);
            else
                zetan=Kn2(k-1)*cot(ve1(k-1))+Kn2(k)*cot(ve2(k))+Kt2(k-1)-Kt2(k);
                zetad=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k)) ...
                    -Kt1(k-1,2)+Kt1(k,1);
            end
            zeta(k)=zetan/zetad;
            zetaaux(k)= zetan/zetad;
        end
    else %Se for um nó do contorno.
        for k=1:nec+1,
            if (k==1)&&(size(P,1)~=size(O,1))
                zetan=Kn2(k)*cot(ve2(k))-Kt2(k);
                zetad=Kn1(k,1)*cot(theta1(k))+Kt1(k,1);
                r(1,1)=1+ (zetan/zetad);
                % comentei porque ja coloquei a norma em Pre_LPEW2 linha 55
                %r(No,1)=(1+ (zetan/zetad))*norm(Qo-T(1,:));
            elseif (k==nec+1)&&(size(P,1)~=size(O,1))
                zetan=Kn2(k-1)*cot(ve1(k-1))+Kt2(k-1);
                zetad=Kn1(k-1,2)*cot(theta2(k-1))-Kt1(k-1,2);
                r(1,2)=1+(zetan/zetad);
                % comentei porque ja coloquei a norma em Pre_LPEW2 linha 55
                %r(No,2)=(1+(zetan/zetad))*norm(Qo-T(nec+1,:));
            else
                zetan=Kn2(k-1)*cot(ve1(k-1))+Kn2(k)*cot(ve2(k))+Kt2(k-1)-Kt2(k);
                zetad=Kn1(k-1,2)*cot(theta2(k-1))+Kn1(k,1)*cot(theta1(k)) ...
                    -Kt1(k-1,2)+Kt1(k,1);
            end
            zeta(k)=zetan/zetad;
        end
    end
    % see eq 25 of the article: "A linearity-preserving cell-centered scheme for the heterogeneous
    %and anisotropic diffusion equations on general meshes" DOI: 10.1002/fld
    for k=1:nec
        if (k==nec)&&(size(P,1)==size(O,1))
            lambda(k)=zeta(k)*Kn1(k,1)*netas(k,1) + zeta(1)*Kn1(k,2)*netas(k,2);

        else
            lambda(k)=zeta(k)*Kn1(k,1)*netas(k,1) + zeta(k+1)*Kn1(k,2)*netas(k,2);

        end
    end

    %% -------------------------------------------------
    % Cálculo de pesos
    %% -------------------------------------------------
    lambda = lambda(:);          % garante vetor coluna
    lambda = lambda(1:nec);      % evita incompatibilidade
    wloc = lambda / sum(lambda); % pesos normalizados
    weight(1,apw(y):apw(y)+nec-1) = wloc;

    apw(y+1) = apw(y) + nec;
    %% -------------------------------------------------
    % Atualiza apw
    %% -------------------------------------------------
    apw(No + 1) = apw(No) + size(O,1);

    %% -------------------------------------------------
    % Interpolação das pressões nos contornos de Neumann
    %% -------------------------------------------------
    if env.config.numcase == 341
        % Vetores locais
        vetor = env.geometry.nsurn1(ns2(No)+1:ns2(No+1));
        comp1 = N(No,1);
        comp2 = N(No,length(vetor));

        % Verifica se o nó pertence ao contorno de Neumann
        if 200 < nflag(No,1) && nflag(No,1) < 300
            % Face comp1
            if env.geometry.bedge(comp1,5) > 200
                a = env.config.bcflag(:,1) == env.geometry.bedge(comp1,5);
                s1 = find(a == 1);
                aa = 0.5*(coord(env.geometry.bedge(comp1,1),:) + coord(env.geometry.bedge(comp1,2),:));
                auxkmap = ferncodes_K(aa(1), aa(2));
                aux1 = r(No,1)*auxkmap(1)*nflagface(s1,2);
            else
                aux1 = 0;
            end

            % Face comp2
            if env.geometry.bedge(comp2,5) > 200
                b = env.config.bcflag(:,1) == env.geometry.bedge(comp2,5);
                s2 = find(b == 1);
                aaa = 0.5*(coord(env.geometry.bedge(comp2,1),:) + coord(env.geometry.bedge(comp2,2),:));
                auxkmap = ferncodes_K(aaa(1), aaa(2));
                aux2 = r(No,2)*auxkmap(1)*nflagface(s2,2);
            else
                aux2 = 0;
            end

            s(No,1) = -(1/sum(lambda))*(aux1 + aux2);
        end
    else
        % Caso numcase ≠ 341
        vetor = env.geometry.nsurn1(ns2(No)+1:ns2(No+1));
        comp1 = N(No,1);
        comp2 = N(No,length(vetor));
        MM = env.geometry.bedge(:,1) == No;
        MMM = find(MM==1);

        if comp1 <= size(env.geometry.bedge,1) && comp2 <= size(env.geometry.bedge,1) && 200 < env.geometry.bedge(MMM,4)
            a = env.config.bcflag(:,1) == env.geometry.bedge(comp1,5);
            s1 = find(a == 1);
            b = env.config.bcflag(:,1) == env.geometry.bedge(comp2,5);
            s2 = find(b == 1);

            s(No,1) = -(1/sum(lambda))*(r(No,1)*env.config.bcflag(s1,2) + ...
                r(No,2)*env.config.bcflag(s2,2));
        end
    end
end

%% -------------------------------------------------
% Saídas
%% -------------------------------------------------
weight = weight(weight ~= 0);
preMPFAD.weight = weight(1,1:apw(nNodes+1)-1);
preMPFAD.s = s;

end

