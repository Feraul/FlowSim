function FS = buildCornerShifts(FS)
%FS.CSR.BUILDCORNERSHIFTS  Precompute per-corner "k / k+1" shift indices.
%
%   FS = fs.csr.buildCornerShifts(FS)
%
%   Given FS.csr populated (cornerNode/cornerElem/nodePtr/nodeNec), add:
%     FS.csr.cornerLocal    [nCorners x 1]  1..nec position within node's corners
%     FS.csr.cornerIsLast   [nCorners x 1]  true when corner is the LAST of its node
%     FS.csr.nodeIsInterior [nNodes   x 1]  true when nns == nec (interior node)
%
%     FS.csr.tCurrent  [nCorners x 1]  flat T (or P) index for T(k) at this corner
%     FS.csr.tNext     [nCorners x 1]  flat T (or P) index for T(k+1) at this corner
%                                       — wraps to T(1) of same node for interior
%                                         last corners; steps to T(k+1) otherwise
%
%   T/P are indexed per NEIGHBOR (length nAllNeighbors = numel(nsurn1)).
%   These per-corner shift arrays let the LPEW2 batched kernels do
%     Tk  = Tall(tCurrent, :)
%     Tk1 = Tall(tNext,    :)
%   without any per-node loop.

    if ~isfield(FS.csr, 'cornerNode') || ~isfield(FS.csr, 'cornerElem')
        error('fs.csr.buildCornerShifts: run fs.csr.buildCorners first');
    end

    nNodes    = FS.mesh.nNodes;
    esurn2    = FS.mesh.esurn2;
    nsurn2    = FS.mesh.nsurn2;

    nec = diff(esurn2);            % [nNodes x 1]  corners per node
    nns = diff(nsurn2);            % [nNodes x 1]  neighbors per node

    % cornerLocal: 1..nec within each node's corner span
    % Efficient vectorization: for each node, (1..nec_i)
    cornerLocal = zeros(FS.csr.nCorners, 1);
    for i = 1:nNodes
        cornerLocal(esurn2(i)+1 : esurn2(i+1)) = 1 : nec(i);
    end

    cornerIsLast   = (cornerLocal == nec(FS.csr.cornerNode));
    nodeIsInterior = (nns == nec);
    isInteriorCorner = nodeIsInterior(FS.csr.cornerNode);

    % T (and P) are indexed per NEIGHBOR — length = numel(nsurn1) = sum(nns)
    % nsurn2 is the row-pointer into T/P for each node's neighbors.
    tBase = nsurn2(FS.csr.cornerNode);       % [nCorners x 1] — 0-based T base for this node

    % tCurrent: flat index of T(k) where k = cornerLocal
    tCurrent = tBase + cornerLocal;          % 1-based

    % tNext: flat index of T(k+1). For interior LAST corner, wrap to T(1) of same node.
    tNext = tCurrent + 1;
    wrapMask = cornerIsLast & isInteriorCorner;
    tNext(wrapMask) = tBase(wrapMask) + 1;   % T(1) of the same node

    % cornerPrev: flat CORNER index (into the corner arrays) of "k-1" within
    % the same node. For interior FIRST corner (k=1), wrap to LAST corner of
    % the same node. For boundary FIRST corner, no wrap — set to sentinel
    % (points to itself, so caller must guard for boundary first corners).
    cornerFlatIdx = (1:FS.csr.nCorners).';
    cornerIsFirst = (cornerLocal == 1);
    cornerBase    = esurn2(FS.csr.cornerNode);   % 0-based corner base for each corner
    cornerPrev = cornerFlatIdx - 1;              % default: previous flat corner
    % For interior first corner: wrap to last corner of same node
    interiorFirstMask = cornerIsFirst & isInteriorCorner;
    cornerPrev(interiorFirstMask) = ...
        cornerBase(interiorFirstMask) + nec(FS.csr.cornerNode(interiorFirstMask));
    % For boundary first corner: sentinel (caller must handle)
    boundaryFirstMask = cornerIsFirst & ~isInteriorCorner;
    cornerPrev(boundaryFirstMask) = cornerFlatIdx(boundaryFirstMask);   % points to self

    % cornerNext: flat CORNER index of "k+1". Interior wraps, boundary steps.
    cornerNext = cornerFlatIdx + 1;
    cornerNext(wrapMask) = cornerBase(wrapMask) + 1;   % first corner of same node
    % For boundary LAST corner: sentinel
    boundaryLastMask = cornerIsLast & ~isInteriorCorner;
    cornerNext(boundaryLastMask) = cornerFlatIdx(boundaryLastMask);

    FS.csr.cornerLocal      = cornerLocal;
    FS.csr.cornerIsFirst    = cornerIsFirst;
    FS.csr.cornerIsLast     = cornerIsLast;
    FS.csr.nodeIsInterior   = nodeIsInterior;
    FS.csr.tCurrent         = tCurrent;
    FS.csr.tNext            = tNext;
    FS.csr.cornerPrev       = cornerPrev;
    FS.csr.cornerNext       = cornerNext;
    FS.csr.boundaryFirstMask = boundaryFirstMask;
    FS.csr.boundaryLastMask  = boundaryLastMask;
end
