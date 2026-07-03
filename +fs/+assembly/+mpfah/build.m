function [M, I, elembedge] = build(env, parms)
%FS.ASSEMBLY.MPFAH.BUILD  Scaffold — delegates to legacy ferncodes_assemblematrixMPFAH.
%
%   PR-C4 scope: creates the +fs/+assembly/+mpfah/ package. Full triplet-form
%   vectorization of the 820-line legacy assembler is a follow-up (PR-C4b).
%   Ships the shape so the OOP contract is uniform across methods.

    pm = env.premethod.MPFAH;
    parameter = pm.parameter;
    weightDMP = pm.weightDMP;
    nflagface = env.config.nflagface;
    dt        = 1;
    SS        = getOr(env.config, 'SS',        0);
    h         = getOr(parms,      'h',         []);
    MM        = getOr(env.config, 'MM',        0);
    gravrate  = getOr(env.preGravity, 'gravrate', 0);
    viscosity = getOr(env.config, 'viscosity', 1);

    [M, I, elembedge] = ferncodes_assemblematrixMPFAH( ...
        parameter, nflagface, weightDMP, SS, dt, h, MM, gravrate, viscosity);
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
