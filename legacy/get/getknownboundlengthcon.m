%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject:  
%Type of file: FUNCTION
%Criate date: 06/07/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%  

%--------------------------------------------------------------------------
%Additional Comments:
%

%--------------------------------------------------------------------------

function [knownboundlength] = getknownboundlengthcon(klb)
%Define global parameters:
global bcflagc normals;

%Initialize "knownboundlength"
knownboundlength = 0;

%It points to a Neumann Boundary Condition 
pointnonnullflag = (bcflagc(:,1) > 200 & bcflagc(:,2) > 0);
flagref = bcflagc(logical(pointnonnullflag),1);
%Choose according "pointnonnullflag"
%There is a Neumann Boundary Condition
if any(klb)
    %Swept the amount of different flags
    for iflag = 1:length(flagref)
        %Swept all edges and get its lengths
        for i = 1:length(klb)
            iedge = klb(i);
            %Attribute the length of each edge
            knownboundlength = knownboundlength + norm(normals(iedge,1:2));
        end  %End of FOR
    end  %End of FOR
%There is NO a Neumann Boundary Condition
else
    knownboundlength = 1;
end  %End of IF
    

