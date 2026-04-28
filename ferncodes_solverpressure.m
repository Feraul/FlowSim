%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in poorus media
%Modified: Fernando Contreras, 2021
function [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq,preMPFAD] = ...
    ferncodes_solverpressure(env,preMPFAD,parmRichardEq,dt,source_wells,tempo)

auxnumcase=env.config.numcase;
auxacel=env.config.acel;
%--------------------------------------------------------------------------
%Solve global algebric system: pressure or hydraulic head
if auxnumcase==331 || (400<auxnumcase && auxnumcase<500)
    %----------------------------------------------------------------------
    % Montagem da matriz global
    [M_old,I_old,] = ferncodes_globalmatrix(env,preMPFAD,parmRichardEq);
    %----------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"

    %Often it may change the global matrix "M" with wells
    [M_old,I_old] = addsource(sparse(M_old),I_old,source_wells,env);

    % Often with source term
    [I_old]=sourceterm(I_old,source_wells);
    %----------------------------------------------------------------------
    % metodo de Picard
    if strcmp(auxacel,'FPI')
        [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq,preMPFAD]=...
       ferncodes_iterpicard(M_old,I_old,preMPFAD,parmRichardEq,env,...
       tempo,dt,source_wells);
        % metodo de Picard com acelaracao de Anderson
    elseif strcmp(auxacel,'AA')
        %% Picard-Anderson Acceleration
        [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq]=...
            ferncodes_iterpicardANLFVPP2(M_old,I_old,preMPFAD,...
            parmRichardEq,env,tempo,dt,source_wells);
    elseif strcmp(auxacel,'Lscheme')
        [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq,preMPFAD]=...
       L_scheme(M_old,I_old,preMPFAD,parmRichardEq,env,...
       tempo,dt,source_wells);
    end
    
else
    % Montagem da matriz global
    [M,I,elembedge] = ferncodes_globalmatrix(env,preMPFAD,parmRichardEq,dt);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    %Often it may change the global matrix "M" with wells
    %[M,I] = addsource(sparse(M),I,source_wells);

    % Often with source term
    [I]=sourceterm(I,source_wells);

    p = solver(M,I);
   
    % auxiliary variables interpolation
    % auxiliary variables interpolation
    [pinterp,~]=ferncodes_pressureinterpNLFVPP(p,preMPFAD,env);
    %Get the flow rate (Diamond)
    [flowrate,flowresult,flowratedif,faceaux] = ferncodes_flowrate(p,pinterp,...
        preMPFAD,parmRichardEq,env);
end





