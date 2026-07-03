function assertFS(FS)
%ASSERTFS  Structural invariant checker for the FS struct.
%
%   assertFS(FS) throws if any of the 5 core invariants fail. Called on every
%   +fs/ solver entry to catch shape drift early.
%
%   Invariants:
%     I1  Top-level fields present: mesh, geom, csr, perm, bc, cfg, state, workspace
%     I2  numel(FS.mesh.esurn2) == FS.mesh.nNodes + 1
%     I3  FS.csr.nCorners == FS.csr.nodePtr(end) - 1  (if csr populated)
%     I4  size(FS.geom.centElem, 1) == FS.mesh.nElems  (if geom populated)
%     I5  size(FS.perm.tensor, 1) == FS.mesh.nElems  (if perm populated)
%
%   Not-yet-populated sub-structs (e.g. empty FS.perm right after mesh.build)
%   are tolerated — invariants only fire when the relevant sub-struct fields exist.

    top = {'mesh','geom','csr','perm','bc','cfg','state','workspace'};
    for k = 1:numel(top)
        if ~isfield(FS, top{k})
            error('assertFS:MissingTopField', 'FS.%s is missing', top{k});
        end
    end

    % I2 — CSR consistency
    if isfield(FS.mesh, 'nNodes') && isfield(FS.mesh, 'esurn2')
        if numel(FS.mesh.esurn2) ~= FS.mesh.nNodes + 1
            error('assertFS:CSRLength', ...
                  'numel(FS.mesh.esurn2)=%d != nNodes+1=%d', ...
                  numel(FS.mesh.esurn2), FS.mesh.nNodes + 1);
        end
    end

    % I3 — corner ptr
    if isfield(FS.csr, 'nCorners') && isfield(FS.csr, 'nodePtr') ...
            && ~isempty(FS.csr.nodePtr)
        if FS.csr.nCorners ~= FS.csr.nodePtr(end) - 1
            error('assertFS:CornerPtr', ...
                  'nCorners=%d != nodePtr(end)-1=%d', ...
                  FS.csr.nCorners, FS.csr.nodePtr(end) - 1);
        end
    end

    % I4 — geom.centElem rows == nElems
    if isfield(FS.geom, 'centElem') && ~isempty(FS.geom.centElem) ...
            && isfield(FS.mesh, 'nElems')
        if size(FS.geom.centElem, 1) ~= FS.mesh.nElems
            error('assertFS:CentElem', ...
                  'size(FS.geom.centElem,1)=%d != nElems=%d', ...
                  size(FS.geom.centElem, 1), FS.mesh.nElems);
        end
    end

    % I5 — perm.tensor rows == nElems
    if isfield(FS.perm, 'tensor') && ~isempty(FS.perm.tensor) ...
            && isfield(FS.mesh, 'nElems')
        if size(FS.perm.tensor, 1) ~= FS.mesh.nElems
            error('assertFS:PermTensor', ...
                  'size(FS.perm.tensor,1)=%d != nElems=%d', ...
                  size(FS.perm.tensor, 1), FS.mesh.nElems);
        end
    end
end
