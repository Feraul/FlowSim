%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media
%Type of file: FUNCTION
%Criate date: 16/04/2015 (My wife is a PHD since yesterday)
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure,flowrate,flowresult,flowratedif] = solvePressure_TPFA(Kde, Kn, ...
    nflagface, Hesq,gravresult,gravrate,gravno,gravelem,fonte,aa,SS,dt,h,MM,wells)
%Define global parameters:

global inedge bedge elem numcase elemarea coord bcflag ...
    gravitational strategy modflowcompared;

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);
I=I+fonte;
% Loop de faces de contorno
m=0;
jj=0;
for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normcont=norm(v0);
    lef=bedge(ifacont,3);
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        % armazena os elementos proximo ao contorno de Dirichlet
         if strcmp(modflowcompared,'y')
            elembedge(jj,1)=bedge(ifacont,3);
            elembedge(jj,2)=nflagface(ifacont,2);
            jj=jj+1;
        end
        c1=nflagface(ifacont,2);
        
        if strcmp(gravitational,'yes')
            if strcmp(strategy,'starnoni')
                m=gravrate(ifacont);
            else
                % proposto de nos
                m1=-nflagface(ifacont,2);
                m=A*(norm(v0)^2*m1-norm(v0)^2*gravelem(lef));
            end
        end
        %Preenchimento
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-c1*A*(norm(v0)^2)+m;
        
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3))- normcont*bcflag(r,2);
        
        
    end
    
end

for iface=1:size(inedge,1),
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(lef, lef)=M(lef, lef)- Kde(iface);
    M(lef, rel)=M(lef, rel)+ Kde(iface);
    M(rel, rel)=M(rel, rel)- Kde(iface);
    M(rel, lef)=M(rel, lef)+ Kde(iface);
    if strcmp(gravitational,'yes')
        if strcmp(strategy,'starnoni')
            m=gravrate(size(bedge,1)+iface);
        else
            m= Kde(iface)*(gravelem(rel,1)-gravelem(lef,1));
            %m=gravrate(size(bedge,1)+iface,1);
            
        end
        I(lef)=I(lef)+m;
        I(rel)=I(rel)-m;
    end
        
end

%--------------------------------------------------------------------
if strcmp(modflowcompared,'y')
    for iw = 1:size(elembedge,1)
        M(elembedge(iw,1),:)=0*M(elembedge(iw,1),:);
        M(elembedge(iw,1),elembedge(iw,1))=1;
        I(elembedge(iw,1))=elembedge(iw,2);
    end
end

% para calcular a carga hidraulica
% para calcular a carga hidraulica
if numcase>300
    if numcase~=336 && numcase~=334 && numcase~=340
        if numcase==333 || numcase==331
            coeficiente=dt^-1*SS.*elemarea(:);
        else
            coeficiente=dt^-1*MM*SS.*elemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            M=M+coeficiente.*eye(size(elem,1));
            I=I+coeficiente.*eye(size(elem,1))*h;
        else
            % Crank-Nicolson method
            I=I+coeficiente.*eye(size(elem,1))*h-0.5*M*h;
            M=0.5*M+coeficiente.*eye(size(elem,1));
            
        end
        
    end
end

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M,I] = addsource(M,I,wells);

%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver".
%This function returns the pressure field with value put in each colocation
%point.
[pressure] = solver(M,I);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%--------------------------------------------------------------------------
%Once the pressure was calculated, the "flowrate" field is also calculated

%Calculate flow rate through edge. "satkey" equal to "1" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[flowrate,flowresult,flowratedif] = calcflowrateTPFA(pressure,Kde,Kn,Hesq,nflagface,1,gravresult,gravrate,aa);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
