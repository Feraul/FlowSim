
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure,flowrate,flowresult,flowratedif] = ...
    ferncodes_solvePressure_TPFA(Kde, Kn, nflag, Hesq,wells,viscosity,...
    Kdec, Knc,nflagc,Con,SS,dt,h,MM,P,time,source)
%Define global parameters:

global inedge bedge elem coord bcflag numcase methodhydro elemarea ...
    modflowcompared Nmod normals;

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);
bedgesize = size(bedge,1);
elembedge=0;
% Loop de faces de contorno
jj=1;
for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(ifacont,1);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        % armazena os elementos proximo ao contorno de Dirichlet
         if strcmp(modflowcompared,'y')
            elembedge(jj,1)=bedge(ifacont,3);
            elembedge(jj,2)=nflag(ifacont,2);
            jj=jj+1;
        end
        %c1=nflag(bedge(ifacont,1),2);
        c1=nflag(ifacont,2);
        
        %Preenchimento
        
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- visonface*A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-visonface*c1*A*(norm(v0)^2);
        
    else
        if numcase==341
            [phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
            %[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
            phi = phiExpNmod10000(1:Nmod);
            C(:,1) = wavenumberExp0Nmod10000(1:Nmod);
            C(:,2) = wavenumberExp1Nmod10000(1:Nmod);
            KMean = 15;
            aaa=0.5*(coord(bedge(ifacont,1),:) + coord(bedge(ifacont,2),:));
            auxkmap = K(aaa(1,1),aaa(1,2),KMean,C(:,1),C(:,2),phi);
            %----------------------------------------------------------
            %auxkmap=kmap(lef, 2);
            I(bedge(ifacont,3)) = I(bedge(ifacont,3))+ normals(ifacont,2)*auxkmap*nflag(ifacont,2);
        else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3))- norm(v0)*bcflag(r,2);
        end
        
    end
    
end

for iface=1:size(inedge,1),
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(bedgesize + iface,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-visonface*Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+visonface*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-visonface*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+visonface*Kde(iface,1);
     
end
if strcmp(modflowcompared,'y')
    for iw = 1:size(elembedge,1)
        M(elembedge(iw,1),:)=0*M(elembedge(iw,1),:);
        M(elembedge(iw,1),elembedge(iw,1))=1;
        I(elembedge(iw,1))=elembedge(iw,2);
    end
end
%--------------------------------------------------------------------------
% para calcular a carga hidraulica
% calcula um problema transiente
if numcase>300
    %
    if numcase~=336 && numcase~=334 && numcase~=335 &&...
            numcase~=337 && numcase~=338 && numcase~=339 &&...
            numcase~=340 && numcase~=341 && numcase~=380 && numcase~=347 && ...
            numcase~=341.1
        if numcase==333 || numcase==331
            %para aquifero nao confinado
            coeficiente=dt^-1*SS.*elemarea(:);
        else
            % para quifero confinado
            coeficiente=dt^-1*MM*SS.*elemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            % equacao 30 Qian et al 2023
            M=coeficiente.*eye(size(elem,1))+M;
            I=I+coeficiente.*eye(size(elem,1))*h;
        else
            % Crank-Nicolson method
            % equacao 33 Qian et al 2023
            I=I+(coeficiente.*eye(size(elem,1))-0.5*M)*h;
            M=  (coeficiente.*eye(size(elem,1))+0.5*M);
        end
        
    end
end
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M,I] = addsource(M,I,wells);

%--------------------------------------------------------------------------
%Solver the algebric system
% Often with source term
[I]=sourceterm(I,source);
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
[flowrate, flowresult,flowratedif]=ferncodes_flowrateTPFA(pressure,Kde,Kn,Hesq,nflag,viscosity,Kdec, Knc,nflagc,Con);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
end
function sol = K(x,y,KMean,C1,C2,phi)
global Nmod varK
coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end
