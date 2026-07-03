function FS = buildCorners(FS)
%BUILDCORNERS  Construct the CSR-flat node→corner layout.
%
%   FS = fs.csr.buildCorners(FS)
%
%   Populates:
%     FS.csr.nCorners      scalar — total corners = sum(esurn2-diff)
%     FS.csr.cornerNode    [nCorners x 1] — which node owns each corner
%     FS.csr.cornerElem    [nCorners x 1] — which element each corner sits in
%     FS.csr.nodePtr       [nNodes+1 x 1] — alias of esurn2 (row-pointer)
%     FS.csr.nodeNec       [nNodes x 1]   — corners per node
%     FS.csr.maxNec        scalar — max corners on any single node
%
%   This layout is the KEY vectorization enabler: with it, every LPEW2 corner
%   computation becomes a flat vector op (see study/artifacts/data-structures.md).

    if ~isfield(FS.mesh, 'esurn1') || ~isfield(FS.mesh, 'esurn2')
        error('fs.csr.buildCorners: FS.mesh.esurn1/esurn2 missing');
    end

    esurn1 = FS.mesh.esurn1(:);
    esurn2 = FS.mesh.esurn2(:);
    nNodes = FS.mesh.nNodes;

    FS.csr.nodePtr    = esurn2;
    FS.csr.nCorners   = numel(esurn1);
    FS.csr.cornerElem = esurn1;
    FS.csr.nodeNec    = diff(esurn2);
    if numel(FS.csr.nodeNec) ~= nNodes
        error('fs.csr.buildCorners: diff(esurn2) has %d entries, expected nNodes=%d', ...
              numel(FS.csr.nodeNec), nNodes);
    end
    FS.csr.cornerNode = repelem((1:nNodes).', FS.csr.nodeNec);
    FS.csr.maxNec     = max(FS.csr.nodeNec);
end
