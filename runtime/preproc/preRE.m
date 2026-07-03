function [parms, env] = preRE(env)
    parms        = env.benchmark.initParms();   % struct vazio padrao
    [parms, env] = env.benchmark.preprocessar(env, parms);
end