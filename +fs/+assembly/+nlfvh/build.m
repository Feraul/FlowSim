function [M, I] = build(env, parms, pinterp)
%FS.ASSEMBLY.NLFVH.BUILD  Scaffold — delegates to legacy ferncodes_assemblematrixNLFVH.
    parameter = env.premethod.MPFAH.parameter;
    viscosity = getOr(env.config, 'viscosity', 1);
    [M, I] = ferncodes_assemblematrixNLFVH(pinterp, parameter, viscosity);
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
