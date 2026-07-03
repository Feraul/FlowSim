%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.  

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [analsw] = hyperb_benchmark(v,t,S0)
%Define global parameters:
global coord centelem elem elemarea numcase;

%Initialize "analsw"
analsw = zeros(size(elem,1),1);

switch numcase
    %Linear Advection obtained from Goosh and Van Altena (2002)
    case 101
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(2.*pi.*(x - v(1).*t)).*sin(2.*pi.*(y - v(2).*t)),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            analsw(i) = (quad2d(@(x,y) ...
                sin(2.*pi.*x).*sin(2.*pi.*y),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2)),...
                'AbsTol',1e-10,'RelTol',1e-10))/elemarea(i);
            
%            analsw(i) = sin(2.*pi*(centelem(i,1) - t)).*sin(2.*pi*(centelem(i,2) - t));
        end  %End of FOR
        
    %Gaussian Hill (Sonar, 2008)
    case 102
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %Define "ref" in each time
            refx = 0.5 + 0.25*cos(t);
            refy = 0.5 + 0.25*sin(t);
            %Obtain the equivalence to "r" (cone)
            r = (quad2d(@(x,y) ((x - refx).^2) + ((y - refy).^2),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            %Attribute the initial condition
            %Define the cylinder with density "1"
            if r <= 0.01
                analsw(i) = 1 - (1/0.01)*r;
            else
                analsw(i) = 0;
            end  %End of IF
        end  %End of FOR
        
    %Sin wave (Wang and Liu, 2004)
    case 103
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(pi.*((x - t) + (y - t))),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            analsw(i) = sin(pi.*((centelem(i,1) - t) + (centelem(i,2) - t)));
        end  %End of FOR
    case 108
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(pi.*((x - t) + (y - t))),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
             analsw(i) = 0.25+0.5*sin(pi.*((centelem(i,1) - t*S0(i,1)) + (centelem(i,2) - t*S0(i,1))));
        end  %End of FOR
    case 107
                
        %Swept all control volumes
        %for i = 1:size(elem,1)
            %Get the vertices for each control volume
            %vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            %vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(pi.*((x - t) + (y - t))),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            %analsw(i) = 0.25+0.5*sin(pi.*((centelem(i,1) - t*S0(i,1)) + (centelem(i,2) - t*S0(i,1))));
            %analsw(i) = (0.5/pi)*sin(2*pi.*((centelem(i,1) - t*S0(i,1))));
            analsw=[0	0;
0.0984772	0.0400957;
0.199915	0.0811225;
0.271837	0.108158;
0.303181	0.118409;
0.336374	0.130528;
0.356656	0.137051;
0.376931	0.141705;
0.393522	0.146362;
0.415642	0.151949;
0.437758	0.156601;
0.465393	0.159379;
0.485647	0.158426;
0.498525	0.154675;
0.501078	-0.156543;
0.515801	-0.159360;
0.530531	-0.160309;
0.545271	-0.158453;
0.563693	-0.156601;
0.580278	-0.153813;
0.598711	-0.149157;
0.618986	-0.144504;
0.642951	-0.137984;
0.670605	-0.129598;
0.698264	-0.120278;
0.740678	-0.104430;
0.783093	-0.0885819;
0.814445	-0.0764617;
0.851331	-0.0615428;
0.901126	-0.0419631;
0.943548	-0.0242457;
1	0];
        %end  %End of FOR
        
    %Two cylinders rotating (Zalezak, 1979)
    case 104
        %Define the angular velocity
        w = [0 0 0.01];
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,:);
            
            %Initialize "velvertex" and "actualvel"
            velvertex = zeros(size(vertcoord,1),3);
            actualvel = zeros(size(vertcoord,1),2);
            
            for j = 1:size(vertcoord,1)
                %Define the velocity for each vertex
                velvertex(j,:) = cross(w,vertcoord(j,:));
                %Calculate the "actualvel"
                actualvel(j,1:2) = vertcoord(j,1:2) - velvertex(j,1:2)*t; 
            end  %End of FOR
            
            %For each element the radianl and angular positions are 
            %calculated
            %Obtain the equivalence to "r1" (cylinder cuted)
            r1 = ((quad2d(@(x,y) (sqrt(x.^2 + (y - 0.25).^2)),...
                min(actualvel(:,1)),max(actualvel(:,1)),...
                min(actualvel(:,2)),max(actualvel(:,2))))/elemarea(i)); 
            %Obtain the equivalence to "r2" (cone)
            r2 = ((quad2d(@(x,y) (sqrt(x.^2 + (y + 0.25).^2)),...
                min(actualvel(:,1)),max(actualvel(:,1)),...
                min(actualvel(:,2)),max(actualvel(:,2))))/elemarea(i)); 

            %Attribute the initial condition
            %Define the cylinder with density "1"
            if r1 < 0.15
                analsw(i) = 1;
            %Define the cone with density variating
            elseif r2 < 0.15
                analsw(i) = 1 - (r2/0.15); 
                
            %In all rest, density receives "0"
            else
                analsw(i) = 0;
            end  %End of IF
    
            %Cut the cylinder
            if (centelem(i,1) < 0.03 && centelem(i,1) > -0.03) && ...
                    (centelem(i,2) > 0.05 && centelem(i,2) < 0.31)
                analsw(i) = 0;
            end  %End of second IF
        end  %End of IF    
    
    %Square walking (Batten et al., 1996)
    case 105
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Define "x" and "y":
            x = centelem(i,1);
            y = centelem(i,2);
            
            %Atribute the value:
            if x > 0.15 + t && x < 0.45 + t && y > 0.15 + t && y < 0.45 + t
                analsw(i) = 1;
            else
                analsw(i) = 0;
            end  %End of IF
        end  %End of FOR
    case 106
        for i = 1:size(elem,1)
            analsw(i) = sin(pi.*((centelem(i,1) - t)));
        end  %End of FOR
end  %End of SWITCH

            