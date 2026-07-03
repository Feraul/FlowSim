
function [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq,preMPFAD]=...
    L_scheme(M_old,RHS_old,...
    preMPFAD,parmRichardEq,env,time,dt,source_wells)

nltol=env.config.nltol;
maxiter=env.config.maxiter;
pmethod=env.config.pmethod;
h_kickoff=parmRichardEq.h_old;
%% Parâmetros do L‑scheme adaptativo
L_min = 1e-3;
L_max = 1e4;
L     = 5;
relax = 1.0;
gamma = 0.9;   % reduz L quando está indo bem
beta  = 2.0;   % aumenta L quando piora

h = parmRichardEq.h_old;          % h^(0)
R_old = norm(M_old*h - RHS_old);  % resíduo inicial
step  = 0;
R0=R_old;
er    = 1;

while ( (er > 1e-7 || maxjump > 1e-3) && step < maxiter )
    step = step + 1;

    n   = size(M_old,1);
    M_L = M_old + L * speye(n);
    RHS_L = RHS_old + L * h;      % usa h^(k) aqui

    % resolve sistema linear
    [p_new, flag] = pcg(M_L, RHS_L, 1e-8, 500);
    if flag ~= 0
        p_new = M_L \ RHS_L;
    end

    h_trial = relax * p_new + (1 - relax) * h;

    dh = h_trial - h;
    maxjump_raw = max(abs(dh));
    maxjump_allowed = 0.2;   % limite físico (0.1–0.5 m)
    scale = min(1, maxjump_allowed / (maxjump_raw + (maxjump_raw==0)));
    h_new = h + scale * dh;
    maxjump = max(abs(h_new - h));
    %fprintf('it=%2d, max|Δh|=%.3e\n', step, maxjump);
    % atualiza estrutura para montar novas matrizes com h_new
    parmRichardEq.h_old = h_new;

    % monta M_new, RHS_new com h_new (seu código)
    
    if strcmp(pmethod,'mpfad')
        [env,parmRichardEq] = PLUG_kfunction(env,parmRichardEq,time);
        [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD);

        % calculo dos pesos que correspondem ao LPEW2
        [preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(preMPFAD,parmRichardEq,env);
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
        % if numcase==439
        %     if 5 <step
        %          dt=1e-4;
        % 
        %      elseif  4<= step && step<= 8
        %         dt=1e-3;
        %      elseif step <8
        %          dt=3e-3;
        %      end
        % parmRichardEq.dt=dt;
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
    R_new = norm(M_new * h_new - RHS_new);

    % erro relativo em relação ao passo anterior
    if R_old > 0
        er_step = R_new / R_old;
    else
        er_step = 0;
    end

    % erro relativo em relação ao inicial (se quiser para parada)
    % erro entre iterações (diagnóstico)
    er_step = R_new / (R_old + (R_old == 0));

    % erro global (critério de parada)
    er = R_new / (R0 + (R0 == 0));

    % ADAPTAÇÃO DE L
     % if R_new > R_old
     %    % resíduo piorou: aumenta L e talvez diminui relax
     %     L = min(L_max, L * beta);
     %     relax = max(0.01, relax * 0.7);
     % else
     %     % resíduo melhorou: reduz L gradualmente
     %     L = max(L_min, L * gamma);
     %     relax = min(1.0, relax * 1.05);
     % end

    % aceita passo
    h    = h_new;
    M_old   = M_new;
    RHS_old = RHS_new;
    R_old   = R_new;
    %fprintf('it=%2d, R=%.3e, er_step=%.3e\n', step, R_new, er_step);
    
    errorelativo(step) = er_step;
end

%--------------------------------------------------------------------------

p=p_new;%M_new\RHS_new;

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