function [p, varargout] = lscheme(varargin)
%FS.ITER.LSCHEME  L-scheme regularized fixed-point iterator (cross-method).
%
%   Wraps legacy L_scheme.m (currently in legacy/unknown/ after PR-F3 move).
%   Called by MetodoTPFA and MetodoMPFAD under case 'LSCHEME'.
    [p, varargout{1:nargout-1}] = L_scheme(varargin{:});
end
