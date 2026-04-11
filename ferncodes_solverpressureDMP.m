
function [p,flowrate,flowresult]=ferncodes_solverpressureDMP(nflagface,...
    parameter,wells,mobility,weightDMP,...
    p_old,w,s,nflagno)
global acel
% interpolação nos nós ou faces
[pinterp]=ferncodes_pressureinterpHP(p_old,nflagface,parameter,weightDMP);
% Montagem da matriz global

[M,I]=ferncodes_assemblematrixDMP(p_old,pinterp,0,parameter,weightDMP,mobility);

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
[p,flowrate,flowresult]=ferncodes_iterpicardANLFVH(M_old,RHS_old,...
    nitpicard,tolpicard,parameter,p_old,nflagface,wells,mobility,weightDMP,w,s,nflagno);
end
end