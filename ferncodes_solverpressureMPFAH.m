function [p,flowrate,flowresult]=ferncodes_solverpressureMPFAH(nflagface,...
             parameter,weightDMP,wells,SS,dt,h,MM,gravrate,viscosity,P,time,source)


[M,I,elembedge]=ferncodes_assemblematrixMPFAH(parameter,nflagface,weightDMP,SS,dt,h,...
    MM,gravrate,viscosity);

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------

% Often with source term
[I]=sourceterm(I,source);
%Solve global algebric system 

% calculo das pressões ou carga hidraulica
p = solver(M,I);

%Message to user:
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpHP(p,nflagface,parameter,weightDMP,0,0,0,0);

%Get the flow rate 
[flowrate,flowresult]=ferncodes_flowratelfvHP(parameter,weightDMP,pinterp,p,viscosity);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
end