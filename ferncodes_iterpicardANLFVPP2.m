function [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq]=...
    ferncodes_iterpicardANLFVPP2(M_old,RHS_old,preMPFAD,...
    parmRichardEq,env,tempo,dt,source_wells)

numcase=env.config.numcase;
pmethod=env.config.pmethod;
elem=env.geometry.elem;

%% calculo do residuo Inicial
if numcase>400
    h_kickoff=parmRichardEq.h_old;
else
    h_kickoff=preNLFV.p_old;
end
R0=norm(M_old*h_kickoff-RHS_old);
zero=zeros(size(env.geometry.elem,1),1);
%% Acelerador de Anderson
%[L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));
[L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-8,'milu','row'));
[p_oldold,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
parmRichardEq.h_old=p_oldold;
%--------------------------------------------------------------------------
[env,parmRichardEq] = PLUG_kfunction(env,parmRichardEq,tempo);
[preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD);
[preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(preMPFAD,parmRichardEq,env);
%--------------------------------------------------------------------------
%[p,erro,iter]=ferncodes_andersonacc2(p_oldold,1e-6,R0,env,parmRichardEq,...
%    preMPFAD,dt,tempo,source_wells);
 [p,erro,iter] = ferncodes_andersonacc2_corrected(p_oldold,1e-6,R0,env,parmRichardEq,...
    preMPFAD,dt,tempo,source_wells);


%Message to user:
fprintf('\n Iteration number, iterations = %d \n',iter)
fprintf('\n Residual error, error = %d \n',erro)
disp('>> The Pressure field was calculated with success!');
%interpolacao das pressoes e concentracoes nos vertices ou faces da malha
%computacional
if strcmp(pmethod,'nlfvpp')
    [pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflagno,w,s,Con,...
        nflagc,wightc,sc);
    %Get the flow rate (Diamond)
    [flowrate,flowresult,flowratedif]=ferncodes_flowrateNLFVPP(p, pinterp,...
        parameter,viscosity,Con,nflagc,wightc,sc,dparameter,cinterp,gravrate);
else
    % auxiliary variables interpolation
    [pinterp,~]=ferncodes_pressureinterpNLFVPP(p,preMPFAD,env);
    %Get the flow rate (Diamond)
    [flowrate,flowresult,flowratedif,faceaux] = ferncodes_flowrate(p,pinterp,...
        preMPFAD,parmRichardEq,env);

end

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end