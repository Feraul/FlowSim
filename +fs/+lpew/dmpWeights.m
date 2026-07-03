function weightDMP = dmpWeights(kmap, elem)
%FS.LPEW.DMPWEIGHTS  DMP-preserving weights (cross-method shared kernel).
%
%   Wraps legacy ferncodes_weightnlfvDMP. Despite the 'nlfv' in the name,
%   this function is called by BOTH MPFA-H and DMP assemblers — it's a
%   general DMP weight computation, not NLFV-specific.
    weightDMP = ferncodes_weightnlfvDMP(kmap, elem);
end
