function [pressure,flowrate,flowresult] = getknownflowrate(elemsize,...
    producelem)
%Define Global parameters:
global bedge inedge normals numcase;

%Initialize "v"
if numcase == 31
    %Initialize "sidelength"
    sidelength = 0.25;
    Q = 1;
    v = [Q 0];
elseif numcase == 31.1
    %Initialize "sidelength"
    sidelength = 75;
    Q = 1.944;%0.02592;
    v = [Q 0];
end  %End of IF

%Initialize "bedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "pressure", "flowrate" and "flowresult"
pressure = ones(elemsize,1);
flowresult = zeros(elemsize,1);
flowrate = zeros(bedgesize + inedgesize,1);

%Swept "bedge"
i = 1:bedgesize;
%Get the flowrate
flowrate(i) = (v(1)*normals(i,1) + v(2)*normals(i,2))/sidelength;
%Attribute "flowrate" to "flowresult"
flowresult(bedge(i,3)) = flowresult(bedge(i,3)) + flowrate(i);

%Swept "inedge"
i = 1:inedgesize;
%Get the flowrate
flowrate(bedgesize + i) = (v(1)*normals(bedgesize + i,1) + ...
    v(2)*normals(bedgesize + i,2))/sidelength;
%Attribute "flowrate" to "flowresult"
flowresult(inedge(i,3)) = ...
    flowresult(inedge(i,3)) + flowrate(bedgesize + i);
flowresult(inedge(i,4)) = ...
    flowresult(inedge(i,4)) - flowrate(bedgesize + i);

%Treats the producer wells
flowresult(producelem) = Q/length(producelem); 