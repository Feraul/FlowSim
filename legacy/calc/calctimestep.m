%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 08/05/2012
%Programer: Fernando R. L. Contreras
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or
%heterogen domain such as isotropic and anisotropic media for each time
%step or in the steady state when will be important.

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

function [dt] = calctimestep(flowrate,satinbound,gamma,Dmedio)
%Define global parameters:
global pormap courant order inedge bedge smethod coord centelem;

%Define the degree of the reconstruction polynomium "n"
n = order - 1;

%Get the flow rate in whole edge when MULTIDIMENSIONAL Schemes are used
%Multidimens. Applic. (Lamine and Edwards, 2010; Kozdon et al.,2011)
if strcmp(smethod,'mwic') || strcmp(smethod,'mwec') || ...
        strcmp(smethod,'rtmd')
    %Join the flowrate calculated for each half-edge in a unic flowrate
    %over the whole edge.
    [flowrate] = joinflowrate(flowrate);
end  %End of IF

%Initialize "dtbyedge"
dtbyedge = zeros(size(inedge,1),1);
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Swept all internal edges
for i = 1:inedgesize
    %Obtain the apropriated deltax:
    %Calculate middle volume: the mean between volume shared by edge
    
    A= norm(coord(inedge(i,1),:)-coord(inedge(i,2),:));
    dx= norm(centelem(inedge(i,3),:)-centelem(inedge(i,4),:));
    
    %Define delta t:
    vel=flowrate(bedgesize+ i)/A;
    dtbyedge(i)= ((courant/((2*n) + 1))*dx^2)/((Dmedio+(abs(vel)*dx)/pormap(1) +gamma*dx^2));
    %dtbyedge(i)= (courant*dx^2)/((Dmedio+(abs(vel)*dx)/pormap(1) +gamma*dx^2));
end  %End of FOR

%--------------------------------------------------------------------------



%Initialize "dtbyboundedge"
dtbyboundedge = zeros(length(satinbound),1);

%Swept edges in "bedge" associated with boundary (injection)
for i = 1:bedgesize
    %Calculate middle volume: the mean between volume shared by edge
  
    lef=bedge(i,3);
    xI = centelem(lef,:);
    xJ=0.5*(coord(bedge(i,1),:)+coord(bedge(i,2),:));
    dx = norm(xI - xJ);
    A=norm(coord(bedge(i,1),:)-coord(bedge(i,2),:)); % comprimento dos baricentros adyacentes a face i
    %Define delta t:
    %Chose according physical effects (gravity existence etc)
    %There is gravity effects
    vel=flowrate(i)/A;
    dtbyboundedge(i)= ((courant/((2*n) + 1))*dx^2)/((Dmedio+(abs(vel)*dx)/pormap(1) +gamma*dx^2));
    %dtbyboundedge(i)= ((courant)*dx^2)/((Dmedio+(abs(vel)*dx)/pormap(1) +gamma*dx^2));
end  %End of FOR

%Do the union between "dtbyedge" and "dtbyboundedge".
dtbyedge = union(dtbyedge,dtbyboundedge);


%Define the values different of "0"
nonzerovalue = logical(dtbyedge ~= 0);
%Finally, define the minor "dt"
dt = min(dtbyedge(nonzerovalue));

