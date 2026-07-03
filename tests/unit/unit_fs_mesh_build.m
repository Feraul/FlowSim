%UNIT_FS_MESH_BUILD  Verify env → FS adapter (fs.mesh.build).
addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_fs_mesh_build');

% Synthesize a minimal env.geometry (skip preprocessormod for isolation)
env.geometry.coord    = [0 0; 1 0; 1 1; 0 1];
env.geometry.elem     = [1 2 3 0 0; 1 3 4 0 0];
env.geometry.bedge    = zeros(4, 5);
env.geometry.inedge   = zeros(1, 6);
env.geometry.nsurn1   = [];
env.geometry.nsurn2   = [0;1;2;3;4];
env.geometry.esurn1   = [1;1;2;1;2;2];
env.geometry.esurn2   = [0;1;3;5;6];
env.geometry.centelem = [1/3 1/3; 2/3 2/3];
env.geometry.elemarea = [0.5; 0.5];
env.geometry.normals  = zeros(4, 3);
env.config.numcase    = 439;
env.config.pmethod    = 'mpfad';
env.config.perm       = [1 0 0 1; 1 0 0 1];
env.config.kmap       = env.config.perm;
env.config.nflag      = zeros(4, 2);
env.config.nflagface  = zeros(4, 2);
env.config.bcflag     = [];

FS = fs.mesh.build(env);

fs_expect(FS.mesh.nNodes  == 4, 'nNodes = 4');
fs_expect(FS.mesh.nElems  == 2, 'nElems = 2');
fs_expect(FS.mesh.nBFaces == 4, 'nBFaces = 4');
fs_expect(FS.mesh.nIFaces == 1, 'nIFaces = 1');
fs_expect(isequal(FS.mesh.coord, env.geometry.coord), 'coord copied');
fs_expect(isequal(FS.geom.centElem, env.geometry.centelem), 'centElem copied');
fs_expect(isequal(FS.geom.elemArea, env.geometry.elemarea), 'elemArea copied');
fs_expect(isfield(FS, 'csr'),       'csr scaffold present');
fs_expect(isfield(FS, 'workspace'), 'workspace scaffold present');
fs_expect(FS.cfg.numcase == 439,    'cfg.numcase copied');
fs_expect(FS.cfg.pmethod == "mpfad" || strcmp(FS.cfg.pmethod, 'mpfad'), 'cfg.pmethod copied');
fs_expect(isequal(FS.perm.tensor, env.config.perm), 'perm.tensor copied');
fs_expect(isequal(FS.bc.nflag,    env.config.nflag), 'bc.nflag copied');

% assertFS on the built FS must pass
try
    fs.util.assertFS(FS);
    fs_expect(true, 'built FS passes assertFS');
catch err
    fs_expect(false, sprintf('built FS failed assertFS: %s', err.message));
end

fs_teardown();
