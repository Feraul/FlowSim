function FS = build(env)
%BUILD  Adapt legacy env (from preprocessormod) → new FS struct.
%
%   FS = fs.mesh.build(env)
%
%   Non-destructive adapter — copies data from env.geometry / env.config
%   into the vectorization-friendly FS layout. Empty sub-structs are
%   created for downstream builders (csr, perm) to fill.
%
%   See ../../my-axon/dev-projects/flowsim-vectorize/phases/study/artifacts/
%       data-structures.md for the full FS spec.

    FS = struct();

    % ── mesh: raw topology ────────────────────────────────────────
    g = env.geometry;
    FS.mesh = struct();
    FS.mesh.nNodes  = size(g.coord,  1);
    FS.mesh.nElems  = size(g.elem,   1);
    FS.mesh.nBFaces = size(g.bedge,  1);
    FS.mesh.nIFaces = size(g.inedge, 1);
    FS.mesh.coord   = g.coord;
    FS.mesh.elem    = g.elem;
    FS.mesh.bedge   = g.bedge;
    FS.mesh.inedge  = g.inedge;
    FS.mesh.nsurn1  = g.nsurn1;
    FS.mesh.nsurn2  = g.nsurn2;
    FS.mesh.esurn1  = g.esurn1;
    FS.mesh.esurn2  = g.esurn2;

    % ── geom: derived, invariant per mesh ─────────────────────────
    FS.geom = struct();
    FS.geom.centElem  = g.centelem;
    FS.geom.elemArea  = g.elemarea;
    if isfield(g, 'normals'), FS.geom.normalInt = g.normals; end

    % ── perm / bc / cfg from env.config (if present) ──────────────
    FS.perm = struct();
    FS.bc   = struct();
    FS.cfg  = struct();
    if isfield(env, 'config')
        c = env.config;
        FS.cfg = c;
        if isfield(c, 'nflag'),     FS.bc.nflag     = c.nflag;     end
        if isfield(c, 'nflagface'), FS.bc.nflagFace = c.nflagface; end
        if isfield(c, 'bcflag'),    FS.bc.bcflag    = c.bcflag;    end
        if isfield(c, 'perm'),      FS.perm.tensor  = c.perm;      end
        if isfield(c, 'kmap'),      FS.perm.kmap    = c.kmap;      end
    end

    % ── csr / state / workspace: empty scaffolds ──────────────────
    FS.csr       = struct();
    FS.state     = struct();
    FS.workspace = struct();
end
