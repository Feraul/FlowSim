function [M, I] = build(env, parms)
%FS.ASSEMBLY.MPFAQL.BUILD  Scaffold — delegates to legacy ferncodes_assemblematrixMPFAQL.
    pm = env.premethod.MPFAQL;
    nflag = env.config.nflag;
    mobility = getOr(env.config, 'mobility', ones(size(env.geometry.inedge, 1), 1));

    [M, I] = ferncodes_assemblematrixMPFAQL( ...
        pm.parameter, pm.weight, pm.s, nflag, pm.weightDMP, mobility);
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
