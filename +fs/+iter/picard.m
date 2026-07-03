function [p, iter, res] = picard(varargin)
%FS.ITER.PICARD  Picard fixed-point iterator (shared across all methods).
%
%   Wraps legacy ferncodes_iterpicard. Called by MetodoTPFA / MetodoMPFAD
%   / MetodoMPFAH / MetodoNLFVPP / MetodoMPFAQL under case 'FPI'.
    [p, iter, res] = ferncodes_iterpicard(varargin{:});
end
