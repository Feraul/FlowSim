tic
[coord,centelem,] = preprocessor;

%Escreve "centelem"
% fid = fopen('C:\Users\Marcio\Doutorado\Programas\Benchmark_Cases\BenchTwophase1_2\centelem64.dat','w'); 
% data1 = centelem(:,1:2)';
% %Print the distribution
% fprintf(fid,'%26.16E %26.16E \r\n',data1);
% fclose(fid);
% pause

%Le a saturação
vecsize = 29498;
%It open the data file ("Start.dat")
readfile = fopen('C:\Users\Marcio\Doutorado\Programas\Benchmark_Cases\BenchTwophase1_2\Saturacao.dat');

getdata = textscan(readfile,'%f\r\n',vecsize);
%Attribute the data to "symaxe"
Sw = cell2mat(getdata);

ponto = 0.64;  %0.72
range = 0.0048;
get_x = centelem(:,1);
get_y = centelem(:,2);

j = 1;
%Swept "centelem"
for i = 1:size(centelem,1)
    if (get_y(i) >= ponto - range) && (get_y(i) <= ponto + range)
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvalue(j) = get_x(i);
        %Attribute to "satfield" the value of "Sw".
        getsatfield(j) = Sw(i);
        %Attribute the number of element to "getelemonline"
        getelemonline(j) = i;
        %Increment "j"
        j = j + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posit = sort(getxvalue);
satfield = zeros(length(getxvalue),1);
elemonline = satfield;
%Reposition the saturation field
for i = 1:length(getxvalue)
    satpointer = logical(getxvalue == posit(i));
    a = getsatfield(satpointer);
    satfield(i) = a(1);
    b = getelemonline(satpointer);
    elemonline(i) = b(1);
end  %End of FOR





%Le a saturação
vecsize = 64*64;
%It open the data file ("Start.dat")
readfile = fopen('C:\Users\Marcio\Doutorado\Programas\Benchmark_Cases\BenchTwophase1_2\centelem64.dat');
getdata = textscan(readfile,'%n %n\r\n',vecsize);
%Attribute the data to "symaxe"
centelem_st = cell2mat(getdata);

%It open the data file ("Start.dat")
readfile3 = fopen('C:\Users\Marcio\Doutorado\Programas\Benchmark_Cases\BenchTwophase1_2\Sw_3aOrdem.dat');

getdata = textscan(readfile3,'%f\r\n',vecsize);
%Attribute the data to "symaxe"
Sw_3 = cell2mat(getdata);

%It open the data file ("Start.dat")
readfile4 = fopen('C:\Users\Marcio\Doutorado\Programas\Benchmark_Cases\BenchTwophase1_2\Sw_4aOrdem.dat');

getdata = textscan(readfile4,'%f\r\n',vecsize);
%Attribute the data to "symaxe"
Sw_4 = cell2mat(getdata);

ponto = 0.64;  %0.72
range = 0.01565;
get_x = centelem_st(:,1);
get_y = centelem_st(:,2);

j = 1;
%Swept "centelem"
for i = 1:size(centelem_st,1)
    if (get_y(i) >= ponto - range) && (get_y(i) <= ponto + range/2)
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvalue3(j) = get_x(i);
        %Attribute to "satfield" the value of "Sw".
        getsatfield3(j) = Sw_3(i);
        %Attribute the number of element to "getelemonline"
        getelemonline3(j) = i;

        %"y_value"
        getxvalue4(j) = get_x(i);
        %Attribute to "satfield" the value of "Sw".
        getsatfield4(j) = Sw_4(i);
        %Attribute the number of element to "getelemonline"
        getelemonline4(j) = i;
        %Increment "j"
        j = j + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posit3 = sort(getxvalue3);
satfield3 = zeros(length(getxvalue3),1);
elemonline3 = satfield3;
%Reposition the saturation field
for i = 1:length(getxvalue3)
    satpointer3 = logical(getxvalue3 == posit3(i));
    a = getsatfield3(satpointer3);
    satfield3(i) = a(1);
    b = getelemonline3(satpointer3);
    elemonline3(i) = b(1);
end  %End of FOR


%Fix "getxvalue" and "getsatfield"
posit4 = sort(getxvalue4);
satfield4 = zeros(length(getxvalue4),1);
elemonline4 = satfield4;
%Reposition the saturation field
for i = 1:length(getxvalue4)
    satpointer4 = logical(getxvalue4 == posit4(i));
    a = getsatfield4(satpointer4);
    satfield4(i) = a(1);
    b = getelemonline4(satpointer4);
    elemonline4(i) = b(1);
end  %End of FOR





plot(posit,satfield,'k','Linewidth',2)
hold on
plot(posit3,satfield3,'b','Linewidth',2)
plot(posit4,satfield4,'r','Linewidth',2)

toc