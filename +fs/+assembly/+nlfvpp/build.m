function [M, I] = build(env, parms, pinterp)
%FS.ASSEMBLY.NLFVPP.BUILD  Scaffold — delegates to legacy ferncodes_assemblematrixNLFVPP.
    pm = env.premethod.NLTPFA;
    parameter = pm.parameter;
    contnorm  = pm.contnorm;
    dt        = 1;
    SS        = getOr(env.config, 'SS',        0);
    h         = getOr(parms,      'h',         []);
    MM        = getOr(env.config, 'MM',        0);
    gravrate  = getOr(env.preGravity, 'gravrate', 0);
    viscosity = getOr(env.config, 'viscosity', 1);
    nflag     = env.config.nflag;

    [M, I] = ferncodes_assemblematrixNLFVPP( ...
        pinterp, parameter, viscosity, contnorm, SS, dt, h, MM, gravrate, nflag);
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
