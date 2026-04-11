%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/11/2015
%Modify data:   /  /2015
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function writes and plots the data calculated.   

%--------------------------------------------------------------------------
%Additional Comments: It is called by the function "validation.m"


%--------------------------------------------------------------------------

function plotandwrite_pressfield(pressure,presanalyt,numcase)

%--------------------------------------------------------------------------
%PLOT RESULTS

%Get the line through the domain and the pressure (cases 1.5 and 1.6)
[pressfield,posit,] = getlineresult(pressure);
%The same for the analytical solution
[pressanal,] = getlineresult(presanalyt);
%Store the results. "posit" is the "r" value and "pressfield" is the
%saturation through this position.
storeresult(posit,pressfield);
    
%Plot the results (Analitical Solution)
figure(3);
plot(posit,pressanal,'LineWidth',2)
%Plot the results (Actual Numerical Solution)
hold on;
plot(posit,pressfield,'-kv');
hold off;

grid on;
title ('Buckley-Leverett Solution')
xlabel ('Reservoir Length');
ylabel ('Water Saturation (pressure)');
    
%Define the limits according to "numcase" number
%Bastian (2002) problem
if numcase == 1.5
    xlim ([0.2 1.2]);
    ylim ([0 1]);
%Any other problem
else
    xlim ([0 5]);
    ylim ([0 8]);
end  %End of IF
%Define a limit for "y" on the graph
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getlineresult"
%--------------------------------------------------------------------------

function [pressfield,posit,elemonline] = getlineresult(pressure)
%Define global parameters:
global centelem;

%Get the number of elements in "z" direction
numelem = round(sqrt(size(centelem,1)));
%Choose a arbitrary position in "centelem"
pointer = ceil(0.5*length(numelem));
%Catch the "y" value in "centelem"
y_value = centelem(numelem(pointer),2);  
%Initialize "j"
j = 1;
%Swept "centelem"
for i = 1:size(centelem,1)
    if (centelem(i,2) >= 0.99*y_value) && ...
            (centelem(i,2) <= 1.01*y_value)
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvalue(j) = centelem(i,1);
        %Attribute to "pressfield" the value of "pressure".
        getsatfield(j) = pressure(i);
        %Attribute the number of element to "getelemonline"
        getelemonline(j) = i;
        %Increment "j"
        j = j + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posit = sort(getxvalue);
pressfield = zeros(length(getxvalue),1);
elemonline = pressfield;
%Reposition the saturation field
for i = 1:length(getxvalue)
    satpointer = logical(getxvalue == posit(i));
    pressfield(i) = getsatfield(satpointer);
    elemonline(i) = getelemonline(satpointer);
end  %End of FOR

%Update "elemonline" without the first positio (it is on face)
elemonline = elemonline(1:length(elemonline));

%--------------------------------------------------------------------------
%Function "storeresult"
%--------------------------------------------------------------------------

function storeresult(pos,pressfield)
%Define global parameters:
global filepath resfolder;

%Create the file name
foldername = resfolder(9:length(resfolder));
filename = [foldername '.dat'];
%Create the file if the scheme is of first order
%Write the file
resfile = fopen([filepath '\' filename],'w'); 
%Print "pos" and "pressfield" values
for i = 1:length(pos)
    fprintf(resfile,'%26.16E %26.16E\r\n',[pos(i)' pressfield(i)]);
end  %End of FOR
%Close the file
fclose(resfile);

