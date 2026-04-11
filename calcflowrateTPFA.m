%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 17/04/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [flowrate,flowresult,flowratedif] = calcflowrateTPFA(p,Kde,Kn,Hesq,nflagface,gravno,gravresult,gravrate,aa)
%Define global parameters
global elem bedge inedge bcflag normals phasekey centelem numcase smethod coord;
flowratedif=0;
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;
mobility=1;
%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
for ifacont=1:size(bedge,1);
    lef=bedge(ifacont,3);
    nor=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*nor);
    if bedge(ifacont,5)<200 % se os nós esteverem na fronteira de DIRICHLET
        c1=nflagface(ifacont,2);
        
        
        flowrate(ifacont)=mobility*A*(nor^2)*(c1-p(lef));
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont)= nor*bcflag(r,2);
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3); %indice do elemento a direita da aresta i
    rel=inedge(iface,4); %indice do elemento a esquerda da aresta i
    
    %-------------------- calculo das vazões e velocidades ---------------%
    
    flowrate(iface+size(bedge,1))=mobility*Kde(iface)*(p(rel)-p(lef));
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
end


%--------------------------------------------------------------------------
%When some multiD schemes are chosen, it is necessary attribute flow rate
%for each half-edge.

%Verify if the scheme is MultiD and which type of that one
if phasekey == 2 && (strcmp(smethod,'mwec') || strcmp(smethod,'mwic') || ...
        strcmp(smethod,'rtmd'))
    %Initialize "auxflowrate"
    auxflowrate = zeros(2*length(flowrate),1);
    %Initialize auxiliary counter
    c = 0;
    %Distribute the flowrate calculated to whole edge in the half-edges.
    for i = 1:length(flowrate)
        auxflowrate(c + 1:c + 2) = 0.5*flowrate(i);
        %Update "c"
        c = c + 2;
    end  %End of FOR
    
    %Finaly, it update "flowrate"
    flowrate = auxflowrate;
    %Clear "auxflowrate"
    clear auxflowrate;
end  %End of IF

