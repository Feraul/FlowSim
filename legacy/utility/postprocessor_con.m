%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 01/05/2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION  
%--------------------------------------------------------------------------

function postprocessor_con(pressure,Concentration,step,normk)
%Define Global parameters:
global coord elem filepath resfolder numcase;

%--------------------------------------------------------------------------
%Create the post file (*.vtk) to PRESSURE Field 
%--------------------------------------------------------------------------
%Initialization parameters

%Attribution of file extension
ext = '.vtk';
%First file name
fname = [filepath '\' resfolder '\'];

%Time level whose the results will appear
step = num2str(step);

%--------------------------------------------------------------------------
%Structure's parameters

%Number of nodes to used
nnode = size(coord,1);
%Number of elements to be used
nelem = size(elem,1);
%"sumnodebyelem" is a parameter which verify how many nodes constitute each
%element multiplied by element amount
sumnodebyelem = sum(sum(elem(:,1:4) ~= 0));
%"countriang" count how many triangle there is in the domain
countriang = sum(elem(:,4) == 0);
%"countquad" count how many quadrangles there is in the domain
countquad = sum(elem(:,4) ~= 0);
%Initialize "type_elem"
type_elem = zeros(size(elem,1),1);

%--------------------------------------------------------------------------
%Write the file

%The command below create the file name joining the parameters which follow:
%"fname" (file name), a tag ("00"), "step" (used in case where is wanted 
%one result by timestep) and the file extension ("step") 
fname_vtk = [fname 'res 00' step ext];
%Open the file to write the text and create the writer "fid".
%The use of letter "w" implies in a quite write, erasing any word in a
%existent file
fid = fopen(fname_vtk,'w'); 
%Write head informations
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data \r\n');
fprintf(fid,'ASCII \r\n');
%Information about grid type
fprintf(fid,'DATASET UNSTRUCTURED_GRID \r\n\r\n');

%Write the POINT informations:
%Head (POINT)
fprintf(fid,'POINTS %i float \r\n\r\n',nnode);
%Distribution (POINT)
data1 = full(coord)';
%Print the distribution
fprintf(fid,'%26.16E %26.16E %26.16E \r\n',data1);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL informations:
%Head (CELL)
fprintf(fid,'CELLS %i %i \r\n\r\n',nelem,sumnodebyelem + nelem);

%There is trianglular elements
if countriang > 0
    %Define how many nodes constitute each element
    %"nodebyelem" is a parameter which verify how many nodes constitute 
    %each element
    nodebyelem = 3;
    %"numdata" is the first column of data to be printed in the CELL sect.
    numdata = nodebyelem*ones(countriang,1);

    %Considering which the initial data (counter) begining from zero, the 
    %node number is diminish of one (-1) in the node definition (command 
    %below). The decision command below is used in order to define if the 
    %elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
    %("nodebyelem" == 4)
    data2 = ...
        [numdata';(elem(1:countriang,1) - 1)';...
        (elem(1:countriang,2) - 1)';(elem(1:countriang,3) - 1)'];
        fprintf(fid,'%i %i %i %i \r\n',data2);
    %Definition of element type (5 is a code to triangle)
    type_elem(1:countriang) = 5;
end  %End of IF (triangles)

%There is quadrangle element
if countquad > 0
    %Define how many nodes constitute each element
    %"nodebyelem" is a parameter which verify how many nodes constitute 
    %each element
    nodebyelem = 4;
    %"numdata" is the first column of data to be printed in the CELL sect.
    numdata = nodebyelem*ones(countquad,1);

    %Considering which the initial data (counter) begining from zero, the 
    %node number is diminish of one (-1) in the node definition (command 
    %below). The decision command below is used in order to define if the 
    %elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
    %("nodebyelem" == 4)
    data2 = ...
        [numdata';(elem(countriang + 1:nelem,1) - 1)';...
        (elem(countriang + 1:nelem,2) - 1)';...
        (elem(countriang + 1:nelem,3) - 1)';...
        (elem(countriang + 1:nelem,4) - 1)'];
    fprintf(fid,'%i %i %i %i %i \r\n',data2);
    %Definition of element type (7 is a code to quadrangle)
    type_elem(countriang + 1:nelem) = 7;            
end  %End of IF (quadrangle)

%Jump a line
fprintf(fid,'\r\n');

%Cell type information
fprintf(fid,'CELL_TYPES %i \r\n\r\n',nelem);
fprintf(fid,'%i \r\n',type_elem);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL_DATA informations:
%Head (CELL_DATA)
fprintf(fid,'CELL_DATA %i \r\n',nelem);

%Write data related to PRESSURE
fprintf(fid,'SCALARS Pressure float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',full(pressure));
%Jump a line
fprintf(fid,'\r\n');

%If the case is a two-phase flow, the oil saturation field is produced
if numcase > 30
    %Write data related to WaterSATURATION
    fprintf(fid,'SCALARS Concentration float 1 \r\n');
    fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
    fprintf(fid,'%26.16E \r\n',Concentration);
    %Jump a line
    fprintf(fid,'\r\n');

end  %End of IF (two-phase flow)

%Write data related to Absolute ERROR
fprintf(fid,'SCALARS NormPermeability float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',normk);
%Jump a line
fprintf(fid,'\r\n');
%Close the text bilder using "fclose"
fclose(fid);

end  %End of function

