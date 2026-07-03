function d = fs_frob(A, B, tol)
%FS_FROB  Relative Frobenius diff between two matrices.
%
%   d = fs_frob(A, B)         % returns ||A-B||_F / ||A||_F
%   d = fs_frob(A, B, 1e-12)  % asserts d < tol via fs_expect
%
%   Handles sparse and dense. Assumes size(A) == size(B).

    if ~isequal(size(A), size(B))
        error('fs_frob: size mismatch (%dx%d vs %dx%d)', ...
              size(A,1), size(A,2), size(B,1), size(B,2));
    end
    denom = norm(A, 'fro');
    if denom == 0
        d = norm(B, 'fro');
    else
        d = norm(A - B, 'fro') / denom;
    end
    if nargin >= 3
        fs_expect(d < tol, sprintf('Frobenius diff %.3e < %.3e', d, tol));
    end
end
