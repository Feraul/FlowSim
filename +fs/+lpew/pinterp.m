function [pressurinterp, cinterp] = pinterp(p, env)
%FS.LPEW.PINTERP  Nodal pressure interpolation (rename of ferncodes_pressureinterpNLFVPP).
%
%   The legacy name is misleading — this function is NOT NLFVPP-specific,
%   it's a general LPEW2-weight-based nodal interpolator used by BOTH
%   MPFA-D and NLFV-PP methods. Moved to +fs/+lpew/ where it belongs.
    [pressurinterp, cinterp] = ferncodes_pressureinterpNLFVPP(p, env);
end
