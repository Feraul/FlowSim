
function [p,flowrate,flowresult,flowratedif,faceaux,parmRichardEq,preMPFAD]=...
    L_scheme(M_old,RHS_old,...
    preMPFAD,parmRichardEq,env,time,dt,source_wells)

nltol=env.config.nltol;
maxiter=env.config.maxiter;
pmethod=env.config.pmethod;
nflag=preMPFAD.nflag;
h_kickoff=parmRichardEq.h_old;
%% Parâmetros do L‑scheme adaptativo
L_min = 1e-3;          % valor mínimo de L
L_max = 1e4;           % valor máximo de L
L     = 10;            % valor inicial (pode vir de uma estimativa)
relax = 1.0;           % fator de relaxação (1.0 = sem relaxação)
gamma = 0.9;           % fator de redução de L quando convergência boa
beta  = 2.0;           % fator de aumento de L quando resíduo sobe

% Estimativa inicial da norma do resíduo (R0 já calculado)
R_old =norm(M_old*h_kickoff-RHS_old);
step=0;
er=1;
while (er > nltol) && (step < maxiter)
    step = step + 1;
    
    % --- L‑scheme com L adaptativo ---
    n = size(M_old,1);
    M_L = M_old + L * speye(n);
    RHS_L = RHS_old + L * parmRichardEq.h_old;
    
    % Resolve o sistema linear (use um solver eficiente, ex: pcg se SPD)
    % p_new = M_L \ RHS_L;           % direto (apenas para pequenos problemas)
    [p_new, flag] = pcg(M_L, RHS_L, 1e-8, 500);  % CG com tolerância
    if flag ~= 0
        warning('PCG não convergiu, usando direto');
        p_new = M_L \ RHS_L;
    end
    
    % Relaxação: h_new = relax * p_new + (1-relax) * h_old
    %h_new = relax * p_new + (1 - relax) * parmRichardEq.h_old;
    
    % Atualiza estrutura (temporariamente)
    parmRichardEq.h_old = p_new;
    
    % --- Atualiza matrizes e RHS (seu código original) ---
    if strcmp(pmethod,'mpfad')
        [env,parmRichardEq] = PLUG_kfunction(env,parmRichardEq,time);
        [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD);
        [preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(preMPFAD,parmRichardEq,env);
        [M_new, I, ~] = ferncodes_globalmatrix(env,preMPFAD,parmRichardEq);
        [M_new, I] = addsource(sparse(M_new), I, source_wells, env);
        [RHS_new] = sourceterm(I, source_wells);
    else
        [pinterp_new,~] = ferncodes_pressureinterpNLFVPP(p_new,nflagno,w,s,Con,...
            nflagc,wightc,sc);
        [M_new, I] = ferncodes_assemblematrixNLFVPP(pinterp_new,parameter,viscosity,...
            contnorm,SS,dt,h,MM,gravrate,nflagno);
        [M_new, RHS_new] = addsource(sparse(M_new), I, wells);
        [RHS_new] = sourceterm(RHS_new, source);
    end
    
    % --- Cálculo do resíduo e adaptação de L ---
    R_new = norm(M_new * p_new - RHS_new);
    if R_old ~= 0
        er = abs(R_new / R_old);
    else
        er = 0;
    end
    
    % Verifica se o resíduo aumentou
   
        % Convergência boa: reduz L gradualmente
        L = max(L_min, L * gamma);
        % Aumenta relaxação de volta para 1
        relax = min(1.0, relax * 1.05);
        % Aceita a solução
        parmRichardEq.h_old = p_new;
        M_old = M_new;
        RHS_old = RHS_new;
        %R_old = R_new;
    
    
    errorelativo(step) = er;
    
    % Opcional: mostrar progresso
    %fprintf('Iter %3d: L=%.2e, res_rel=%.2e, relax=%.2f\n', step, L, er, relax);
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