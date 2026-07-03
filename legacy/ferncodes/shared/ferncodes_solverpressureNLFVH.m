function [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_solverpressureNLFVH(nflagface,...
    parameter,wells,mobility,weightDMP,...
    p_old,w,s,nflagno,contnorm,weightDMPc,Con,nflagfacec,dparameter,wightc,sc,SS,dt,h,MM,gravrate)


global acel

% interpolação nos nós ou faces
[pinterp,]=ferncodes_pressureinterpHP(p_old,nflagface,parameter,weightDMP,...
                      weightDMPc,Con,nflagfacec,dparameter);

[M,I]=ferncodes_assemblematrixNLFVH(pinterp,parameter,mobility);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);
% Picard iteration
if strcmp(acel,'FPI')
    [p,flowrate,flowresult]=ferncodes_iterpicardNLFVH(M_old,RHS_old,...
        nitpicard,tolpicard,parameter,p_old,nflagface,wells,mobility,weightDMP);
    % Iteração de Picard com acelerador de Anderson
elseif strcmp(acel,'AA')
    [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_iterpicardANLFVH(M_old,RHS_old,...
        parameter,p_old,nflagface,wells,mobility,weightDMP,w,s,...
        nflagno,contnorm,weightDMPc,Con,nflagfacec,dparameter,wightc,sc,SS,dt,h,MM,gravrate);
end
end