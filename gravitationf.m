function [g]=gravitationf(kmap,gravelem)
global inedge bedge elem centelem coord dens
Klef=zeros(2,2);
Krel=zeros(2,2);
G=zeros(size(elem,1),1);
R=[0 1 ;-1 0 ];
K1=zeros(3,3);
K2=zeros(3,3);
K=zeros(3,3);
for ifacont=1:size(bedge,1)
    % elemento a esquerda 
    lef=bedge(ifacont,3);
    % vetor de orientacao da face em questao
    ve1=coord(bedge(ifacont,2),1:2)-coord(bedge(ifacont,1),1:2);
    
    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
    
    Keq=Klef;
    g(ifacont,1)=dens(1,1)*dot((R*ve1')'*Keq,(gravelem(lef,1:2)));
    G(lef,1)=G(lef,1)-g(ifacont,1);

end
for iface=1:size(inedge,1)
         % elementos a esquerda e a direita
         lef=inedge(iface,3);
         rel=inedge(iface,4);
         
         % calculo do ponto meio da face
         vm=(coord(inedge(iface,1),1:2)+coord(inedge(iface,2),1:2))*0.5;
        % vd1aux=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
         vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
         % calculo da distancia do centro ao ponto meio da face
         dj1=norm(centelem(lef,1:2)-vm);
         dj2=norm(centelem(rel,1:2)-vm);

         % tensor de permeabilidade do elemento a esquerda
         
         Klef(1,1)=kmap(elem(lef,5),2);
         Klef(1,2)=kmap(elem(lef,5),3);
         Klef(2,1)=kmap(elem(lef,5),4);
         Klef(2,2)=kmap(elem(lef,5),5);
         
         % tensor de permeabilidade do elemento a direita
         
         Krel(1,1)=kmap(elem(rel,5),2);
         Krel(1,2)=kmap(elem(rel,5),3);
         Krel(2,1)=kmap(elem(rel,5),4);
         Krel(2,2)=kmap(elem(rel,5),5);
         
         Keq=inv((dj1*inv(Klef)+dj2*inv(Krel))); % equation 21
         graveq=((dj1*gravelem(lef,:)+dj2*gravelem(rel,:))'); % equation 22
         g(iface+size(bedge,1),1)=dens(1,1)*dot(((R*vd1')')*Keq, graveq);% equation 20
        
         G(lef,1)=G(lef,1)-g(iface+size(bedge,1),1);
         G(rel,1)=G(rel,1)+g(iface+size(bedge,1),1);
     
end   
end
