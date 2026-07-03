function [M, I] = build(env, parms, p, pinterp, gamma, graviterm)
%FS.ASSEMBLY.DMP.BUILD  Scaffold — delegates to legacy ferncodes_assemblematrixDMP.
    pm = env.premethod.MPFAH;
    parameter = pm.parameter;
    weightDMP = pm.weightDMP;
    mobility  = getOr(env.config, 'mobility', ones(size(env.geometry.inedge, 1), 1));
    if nargin < 5, gamma = []; end
    if nargin < 6, graviterm = []; end
    [M, I] = ferncodes_assemblematrixDMP( ...
        p, pinterp, gamma, parameter, weightDMP, mobility, graviterm);
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
