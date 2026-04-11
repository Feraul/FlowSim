function [gravrate,gravresult]=gravitation(kmap)
global inedge bedge elem centelem coord
K1=zeros(3,3);
K2=zeros(3,3);
K=zeros(3,3);
gravresult=zeros(size(elem,1),1);
R=[0 1 ;-1 0 ];

for ifacont=1:size(bedge,1)
    % elemento a esquerda
    lef=bedge(ifacont,3);
    % vetor de orientacao da face em questao
    ve1=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    vm=0.5*(coord(bedge(ifacont,2),:)+coord(bedge(ifacont,1),:));
    % tensor de permeabilidade do elemento a esquerda
    % Keq(1,1)=kmap(elem(lef,5),2);
    % Keq(1,2)=kmap(elem(lef,5),3);
    % Keq(2,1)=kmap(elem(lef,5),4);
    % Keq(2,2)=kmap(elem(lef,5),5);
    C1 = centelem(bedge(ifacont,3),:);
    %Determinação das alturas dos elementos à esquerda.
    % ve1 = coord(bedge(ifacont,2),:) - coord(bedge(ifacont,1),:); % face
    ve2 = coord(bedge(ifacont,2),:) - C1; %Do centro esquerdo ao fim da face.

    ce = cross(ve1,ve2); % produto vetorial
    Hesq1 = norm(ce)/norm(ve1); % altura a relativo as faces do contorno


    %Essa é UMA maneira de construir os tensores
    K(1,1) = kmap(elem(bedge(ifacont,3),5),2);
    K(1,2) = kmap(elem(bedge(ifacont,3),5),3);
    K(2,1) = kmap(elem(bedge(ifacont,3),5),4);
    K(2,2) = kmap(elem(bedge(ifacont,3),5),5);

    %Cálculo das constantes tangenciais e normais

    Kn1 = (RotH(ve1)'*K*RotH(ve1))/norm(ve1)^2;
    A=-Kn1/(Hesq1);
    %Keq=Klef;
    if bedge(ifacont,5)>200

        %gravrate(ifacont,1)=dot((R*ve1')'*Keq,([0 1 0]));
        gravrate(ifacont,1)=0;
        %gravrate(ifacont,1)=A*(norm(ve1))*(vm(1,2)-centelem(lef,2));
    else
        % if bedge(ifacont,5)==102
        %   gravrate(ifacont,1)=dot((R*ve1')'*Keq,([0 centelem(lef,2)]));
        % else
        % gravrate(ifacont,1)=dot((R*ve1')'*Keq,([0 1 0]));
        % end


        gravrate(ifacont,1)=A*(norm(ve1))*(vm(1,2)-centelem(lef,2));
    end
    gravresult(lef,1)=gravresult(lef,1)-gravrate(ifacont,1);

end
for iface=1:size(inedge,1)
    % elementos a esquerda e a direita
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %
    %   % calculo do ponto meio da face
    %   vm=(coord(inedge(iface,1),1:2)+coord(inedge(iface,2),1:2))*0.5;
    %  % vd1aux=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    %   vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
    %   % calculo da distancia do centro ao ponto meio da face
    %   dj1=norm(centelem(lef,1:2)-vm);
    %   dj2=norm(centelem(rel,1:2)-vm);
    %
    %   % tensor de permeabilidade do elemento a esquerda
    %
    %  R1=[0 1 ;-1 0 ];
    %  K3(1,1)=kmap(elem(lef,5),2);
    %  K3(1,2)=kmap(elem(lef,5),3);
    %  K3(2,1)=kmap(elem(lef,5),4);
    %  K3(2,2)=kmap(elem(lef,5),5);
    %
    %  % tensor de permeabilidade do elemento a direita
    %  K4(1,1)=kmap(elem(rel,5),2);
    %  K4(1,2)=kmap(elem(rel,5),3);
    %  K4(2,1)=kmap(elem(rel,5),4);
    %  K4(2,2)=kmap(elem(rel,5),5);
    %
    %  Keq=inv((dj1*inv(K3)+dj2*inv(K4))); % equation 21
    %  graveq=((dj1*[0 1]+dj2*[0 1])'); % equation 22
    % % gravrate(iface+size(bedge,1),1)=dot(((R1*vd1')'),Keq*graveq);% equation 20
    %
    %  gravrate(iface+size(bedge,1),1)=norm(vd1)*Keq*(centelem(lef,2)-centelem(rel,2));% equation 20
    %

    %Determinação dos centróides dos elementos à direita e à esquerda.%
    C1 = centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
    C2 = centelem(inedge(iface,4),:); % baricentro do elemento direito
    vcen = C2 - C1;
    vd1 = coord(inedge(iface,2),:) - coord(inedge(iface,1),:);

    %Determinação das alturas dos centróides dos elementos à direita e à%
    %esquerda.                                                          %

    vd2 = C2 - coord(inedge(iface,1),:);     %Do início da aresta até o
    %centro da célula da direita.
    cd = cross(vd1,vd2);
    H2 = norm(cd)/norm(vd1); % altura a direita

    ve2 = C1 - coord(inedge(iface,1),:);

    ce = cross(vd1,ve2);
    H1 = norm(ce)/norm(vd1); % altura a esquerda

    %Cálculo das constantes.%
    %A segunda entrada será tal que: 1=dir, 2=esq.

    %Essa é UMA maneira de construir os tensores.
    %Permeability on the Left

    K1(1,1) = kmap(elem(inedge(iface,3),5),2);
    K1(1,2) = kmap(elem(inedge(iface,3),5),3);
    K1(2,1) = kmap(elem(inedge(iface,3),5),4);
    K1(2,2) = kmap(elem(inedge(iface,3),5),5);

    %Permeability on the Right

    K2(1,1) = kmap(elem(inedge(iface,4),5),2);
    K2(1,2) = kmap(elem(inedge(iface,4),5),3);
    K2(2,1) = kmap(elem(inedge(iface,4),5),4);
    K2(2,2) = kmap(elem(inedge(iface,4),5),5);

    % calculo das constantes tangenciais e normais em cada face interna
    Kn1 = (RotH(vd1)'*K1*RotH(vd1))/norm(vd1)^2;
    Kn2 = (RotH(vd1)'*K2*RotH(vd1))/norm(vd1)^2;
    % calculo das constantes nas faces internas
    Kde = -norm(vd1)*((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
    gravrate(iface+size(bedge,1),1)=Kde*(centelem(rel,2)-centelem(lef,2));% equation 20


    gravresult(lef,1)=gravresult(lef,1)-gravrate(iface+size(bedge,1),1);
    gravresult(rel,1)=gravresult(rel,1)+gravrate(iface+size(bedge,1),1);

end
end
