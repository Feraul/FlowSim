function [ O, P, T, Qo ] = OPT_Interp_LPEW (No,env)

%Retorna os vetores O, P, T e Qo.
% Lembrando que estes esurn1, nsurn1 já estan ordenados em sentido
% anti-horario, sequencialmente. 

%Pré-alocação dos vetores.%

P=zeros(env.geometry.nsurn2(No+1)-env.geometry.nsurn2(No),3); % vetor de pontos na vizinhança do nó "ni".
T=zeros(env.geometry.nsurn2(No+1)-env.geometry.nsurn2(No),3); % vetor de pontos dinamicos na vizinhança do nó "ni".
O=zeros(env.geometry.esurn2(No+1)-env.geometry.esurn2(No),3); % vetor de baricentro na vizinhança do nó "ni".
Qo=env.geometry.coord(No,:);                     % coordenada do nó "ni".



%Construção dos vetores P, dos nós vizinhos ao nó "ni", e T, dos pontos%
%médios das fases que concorrem no nó "ni".                            %
for i=1:size(P,1),
    P(i,:)=env.geometry.coord(env.geometry.nsurn1(env.geometry.nsurn2(No)+i),:);
    T(i,:)=(P(i,:)+Qo)/2;
end

%Construção do vetor O, dos centróides (pontos de colocação) dos elementos%
%que concorrem no nó ni.                                                  %

for i=1:size(O,1)
    %Verifica se o elemento é um quadrilátero ou um triângulo.
    if env.geometry.elem(env.geometry.esurn1(env.geometry.esurn2(No)+i),4)==0 % lenbrando que o quarta columna
        b=3;                  
    else
        b=4;  % da matriz de elementos é para quadrilateros
    end
    %Carrega adequadamente o vetor O (braicentro de cada elemento)
    for j=1:b
        O(i,1)=O(i,1)+(env.geometry.coord(env.geometry.elem(env.geometry.esurn1(env.geometry.esurn2(No)+i),j),1)/b);
        O(i,2)=O(i,2)+(env.geometry.coord(env.geometry.elem(env.geometry.esurn1(env.geometry.esurn2(No)+i),j),2)/b);
        O(i,3)=O(i,3)+(env.geometry.coord(env.geometry.elem(env.geometry.esurn1(env.geometry.esurn2(No)+i),j),3)/b);
    end
end

end
