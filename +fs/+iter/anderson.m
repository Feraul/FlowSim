function [x, iter, res_hist] = anderson(varargin)
%FS.ITER.ANDERSON  Anderson-accelerated fixed-point iterator (cross-method).
%
%   Wraps legacy ferncodes_andersonacc. Historically its file lived in
%   legacy/ferncodes/shared/ but the name suggests NLFVPP-specificity —
%   it's actually cross-method (used by MPFA-D 'AA' path too). This
%   thin wrapper documents the true home; callers using +fs/ invoke
%   fs.iter.anderson; legacy callers keep using the ferncodes_ name.
    [x, iter, res_hist] = ferncodes_andersonacc(varargin{:});
end
