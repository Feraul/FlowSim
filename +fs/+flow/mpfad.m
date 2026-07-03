function [flowrate, flowresult, flowratedif, faceaux] = mpfad(p, pinterp, env)
%FS.FLOW.MPFAD  Vectorized MPFA-D flow-rate (delegates to legacy ferncodes_flowrate).
%
%   Legacy ferncodes_flowrate.m is already vectorized (marked "VETORIZADO"
%   in header). This is a rename/move into the +fs/+flow/ package.
    [flowrate, flowresult, flowratedif, faceaux] = ferncodes_flowrate(p, pinterp, env);
end
