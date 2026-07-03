
%Programer: Márcio Souza
%Modified: Fernando Contreras,
function [p,flowrate,flowresult,flowratedif,flowresultc]=...
                                ferncodes_solverpressureNLFVPP(nflag,...
                                  parameter,kmap,wells,viscosity,V,N,...
                                  p_old,contnorm,weight,s,Con,nflagc,...
                                  weightc,sc,weightDMPc,nflagfacec,dparameter,...
                                  SS,dt,h,MM,gravrate,source)
                            
%Define global parameters
global acel;

% interpolação das pressoes e concentracoes nos vetrtices ou faces
[pinterp,]=ferncodes_pressureinterpNLFVPP(p_old,nflag,weight,s,Con,nflagc,...
                                                               weightc,sc);

% montagem das matriz global 
[M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,contnorm,...
                                                SS,dt,h,MM,gravrate,nflag);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);
% Often with source term
[RHS_old]=sourceterm(RHS_old,source);
%% full Picard iteration
if strcmp(acel,'FPI')
    [p,flowrate,flowresult,flowratedif]=ferncodes_iterpicard(M_old,RHS_old,...
        parameter,weight,s,p_old,nflag,wells,viscosity,Con,nflagc,weightc,sc,...
        dparameter,contnorm,SS,dt,h,MM,gravrate,source);
elseif strcmp(acel,'AA')
    %% Picard-Anderson Acceleration
    [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
        parameter,weight,s,p_old,nflag,wells,viscosity,0,contnorm,Con,nflagc,...
        weightc,sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate,source);
end
end