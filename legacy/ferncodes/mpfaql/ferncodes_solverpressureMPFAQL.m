
%Modified: Fernando Contreras, 2021

function [p,flowrate,flowresult]=ferncodes_solverpressureMPFAQL(nflagno,...
    parameter,kmap,weightDMP,wells,mobility,V,Sw,N,weight,s)

[M,I]=ferncodes_assemblematrixMPFAQL(parameter,weight,s,nflagno,weightDMP,mobility);

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p = solver(M,I);

%Message to user:
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpMPFAQL(p,nflagno,weight,s);
%Get the flow rate (Diamond)
[flowrate,flowresult]=ferncodes_flowratelfvMPFAQL(parameter,weightDMP,mobility,pinterp,p);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end