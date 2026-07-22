%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in poorus media
%Modified: Fernando Contreras, 2021

function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
    ferncodes_solver(env, parms, dt, source_wells, tempo)

    % ── Monta matriz global — delega ao metodo ────────────────────
    [M, I] = env.metodo.montarSistema(env, parms, dt);

    % ── Adiciona fontes ───────────────────────────────────────────
    [M, I] = addsource(sparse(M), I, source_wells, env); %poços
    [I]    = sourceterm(I, source_wells);

    % ── Resolve — delega ao metodo (LINEAR, Picard, AA, L-scheme) ────────
    [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
        env.metodo.resolver(M, I, parms, env, tempo, dt, source_wells);
end
