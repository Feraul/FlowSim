function [p,flowrate,flowresult,flowratedif,faceaux,parms,env] = ...
        ferncodes_iterpicard(M_old, RHS_old, parms, env, time, dt, source_wells)

    nltol   = env.config.nltol;
    maxiter = env.config.maxiter;

    p_old = parms.h_old;
    R0    = norm(M_old*p_old - RHS_old);

    step = 0;
    er   = 1;

    while (er >= nltol) && (step < maxiter)
        step = step + 1;

        p_new       = solver(M_old, RHS_old);
        p_new       = p_old + 0.5*(p_new - p_old);
        parms.h_old = p_new;

        % ── atualiza K(h) SOMENTE se o modelo fisico exigir ──────
        % (Richards: sim, precisa. NL-TPFA/MPFA com K constante: nao precisa,
        %  a nao-linearidade esta so no esquema numerico, nao na fisica)
        if env.benchmark.precisaAtualizarPermeabilidade()
            [env, parms] = PLUG_kfunction(env, parms, time);
            [env]        = env.metodo.atualizarPremethod(env, parms);
        end

        % ── remonta o sistema (sempre necessario -- mesmo com K fixo,
        %    a matriz pode depender nao-linearmente de p via o esquema) ──
        [M, I] = env.metodo.montarSistema(env, parms, dt);

        [M_new, I] = addsource(sparse(M), I, source_wells, env);
        [RHS_new]  = sourceterm(I, source_wells);

        R = norm(M_new*p_new - RHS_new);
        if R0 ~= 0.0
            er = abs(R/R0);
        else
            er = 0.0;
        end
        errorelativo(step) = er;

        M_old   = M_new;
        RHS_old = RHS_new;
        p_old   = p_new;
    end

    p = p_new;

    fprintf('\n Iteration number, iterations = %d \n', step);
    fprintf('\n Residual error, error = %d \n', er);
    disp('>> The Pressure field was calculated with success!');

    [flowrate, flowresult, flowratedif, faceaux] = ...
        env.metodo.calcularFlowrate(p, env, parms);

    disp('>> The Flow Rate field was calculated with success!');
end