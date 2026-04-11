function [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq]=...
    ferncodes_iterpicard(M_old,RHS_old,...
    preMPFAD,parmRichardEq,env,time,dt,source_wells)

% incialicando parametros globais
nltol=env.config.nltol;
maxiter=env.config.maxiter;
pmethod=env.config.pmethod;
nflag=preMPFAD.nflag;
h_kickoff=parmRichardEq.h_old;
%% calculo do residuo Inicial
R0=norm(M_old*h_kickoff-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;
zero=zeros(size(env.geometry.elem,1),1);
while (nltol<er || nltol==er) && (step<maxiter)
    % atualiza iterações
    step=step+1;

    % calculo das pressões
    % [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-8));

    %  [p_new,]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);

    p_new = solver(M_old,RHS_old);
    parmRichardEq.h_old=p_new;
    if strcmp(pmethod,'mpfad')
        [env,parmRichardEq] = PLUG_kfunction(env,parmRichardEq,time);
        [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD,time);
        
        % calculo dos pesos que correspondem ao LPEW2
        [preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(zero,preMPFAD,parmRichardEq,env);
        %==================================================================
        % if numcase==432
        %     if 7 <step
        %         dt=0.7*dt;
        %     elseif  3<= step && step<= 7
        %         dt=1*dt;
        %     elseif step <7
        %         dt=1.3*dt;
        %     end
        % end
        %==================================================================
        % Montagem da matriz global
        [M,I,] = ferncodes_globalmatrix(env,preMPFAD,parmRichardEq);
        %------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        %Often it may change the global matrix "M" with wells
        [M_new,I] = addsource(sparse(M),I,source_wells,env);

        % Often with source term
        [RHS_new]=sourceterm(I,source_wells);
    else
        %% plotagem no visit
        %S=ones(size(p_new,1),1);
        %ferncodes_postprocessor(p_new,S,step)
        [pinterp_new,]=ferncodes_pressureinterpNLFVPP(p_new,nflagno,w,s,Con,...
            nflagc,wightc,sc);
        %% Calculo da matriz global

        [M,I]=ferncodes_assemblematrixNLFVPP(pinterp_new,parameter,viscosity,...
            contnorm,SS,dt,h,MM,gravrate,nflagno);
        %--------------------------------------------------------------------------
        %Often it may change the global matrix "M"
        [M_new,RHS_new] = addsource(sparse(M),I,wells);
        % Often with source term
        [RHS_new]=sourceterm(RHS_new,source);
    end
    %% Calculo do residuo

    R = norm(M_new*p_new - RHS_new);

    if (R0 ~= 0.0)
        er = abs(R/R0);
    else
        er = 0.0; %exact
    end
    errorelativo(step)=er;

    % atualizar
    M_old=M_new;
    RHS_old=RHS_new;
    %h=p_new;

end
dt_aux=dt;
%--------------------------------------------------------------------------
p=p_new;

%Message to user:
fprintf('\n Iteration number, iterations = %d \n',step)
fprintf('\n Residual error, error = %d \n',er)
%Message to user:
disp('>> The Pressure field was calculated with success!');
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