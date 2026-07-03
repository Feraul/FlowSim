%UNIT_ASSERTFS  Verify FS invariant checker fires on the right shapes.
addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_assertFS');

% ── Build a minimal valid FS ──────────────────────────────────────
FS = struct();
FS.mesh = struct('nNodes', 4, 'nElems', 2, 'nBFaces', 4, 'nIFaces', 1, ...
                 'coord', [0 0; 1 0; 1 1; 0 1], 'elem', [1 2 3 0; 1 3 4 0], ...
                 'bedge', zeros(4,5), 'inedge', zeros(1,6), ...
                 'nsurn1', [], 'nsurn2', [0;1;2;3;4], ...
                 'esurn1', [1;1;2;1;2;2], 'esurn2', [0;1;3;5;6]);
FS.geom      = struct('centElem', [1/3 1/3; 2/3 2/3], 'elemArea', [0.5; 0.5]);
FS.csr       = struct();  % empty ok
FS.perm      = struct();  % empty ok
FS.bc        = struct();
FS.cfg       = struct();
FS.state     = struct();
FS.workspace = struct();

% Should PASS
try
    fs.util.assertFS(FS);
    fs_expect(true, 'valid FS passes assertFS');
catch err
    fs_expect(false, sprintf('valid FS unexpectedly failed: %s', err.message));
end

% ── Break I1 (missing top field) ──────────────────────────────────
bad = rmfield(FS, 'perm');
gotErr = false;
try, fs.util.assertFS(bad); catch, gotErr = true; end
fs_expect(gotErr, 'missing top field triggers assertion');

% ── Break I2 (esurn2 length mismatch) ─────────────────────────────
bad = FS;
bad.mesh.esurn2 = [0;1;3;5];  % length 4 but nNodes=4 → needs 5
gotErr = false;
try, fs.util.assertFS(bad); catch, gotErr = true; end
fs_expect(gotErr, 'esurn2 length mismatch triggers assertion');

% ── Break I4 (centElem row count mismatch) ────────────────────────
bad = FS;
bad.geom.centElem = [1/3 1/3];  % 1 row but nElems=2
gotErr = false;
try, fs.util.assertFS(bad); catch, gotErr = true; end
fs_expect(gotErr, 'centElem row-count mismatch triggers assertion');

% ── Break I5 (perm.tensor row count mismatch) ─────────────────────
bad = FS;
bad.perm.tensor = [1 0 0 1];  % 1 row but nElems=2
gotErr = false;
try, fs.util.assertFS(bad); catch, gotErr = true; end
fs_expect(gotErr, 'perm.tensor row-count mismatch triggers assertion');

fs_teardown();
