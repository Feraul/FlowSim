function [flowrate, flowresult, flowratedif, faceaux] = tpfa(p, env)
%FS.FLOW.TPFA  Vectorized TPFA flow-rate (delegates to legacy ferncodes_flowrateTPFA).
    [flowrate, flowresult, flowratedif, faceaux] = ferncodes_flowrateTPFA(p, env);
end
