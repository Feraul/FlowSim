%It is called by "preMPFA.m"

function [env] = ferncodes_calflag(env,parmRichardEq,time)
  

    nelem_nodes = size(env.geometry.coord,  1);
    nelem_faces = size(env.geometry.bedge,  1);

    % Inicializa
    nflag     = 5000 * ones(nelem_nodes, 2);
    nflagface = zeros(nelem_faces, 2);

    % Injeta no env para o benchmark usar
    env.config.nflag_init     = nflag;
    env.config.nflagface_init = nflagface;

    % ── Delega ao benchmark ───────────────────────────────────────
    [nflag, nflagface] = env.benchmark.configurarFlags(env, parmRichardEq, time);

    env.config.nflag     = nflag;
    env.config.nflagface = nflagface;
end
