filepath = 'C:\Users\Marcio\Doutorado\Programas\Benchmark_MeshFiles';
filename_entrada = 'pontosinternos.dat';
filename_saida = 'pontosinternos2.dat';
%It open the data file ("Start.dat")
readfile = fopen([filepath '\' filename_entrada]);
getdata = textscan(readfile,'%f %f %f %f',2500);
num = cell2mat(getdata);

a = -1;
b = 1;
rx = a + (b-a).*rand(2500,1);
ry = a + (b-a).*rand(2500,1);


num(:,2) = num(:,2) + rx*0.32*(1/51); 
num(:,3) = num(:,3) + ry*0.32*(1/51); 

resfile = fopen([filepath '\' filename_saida],'w'); 
%Print "pos" and "satfield" values
for i = 1:2500
    fprintf(resfile,'%d %26.16E %26.16E %d\r\n',num(i,:));
end  %End of FOR
%Close the file
fclose(resfile);
