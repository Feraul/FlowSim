addpath(fullfile(pwd, 'tests', 'helpers'));  %% bootstrap so fs_* helpers resolve
%SMOKE_MESH_LOAD  Load every .msh fixture and report node + element counts.
%
%   Verifies mesh files are readable and roughly sized as expected.
%   Uses a lightweight parser — we count `$Nodes` / `$Elements` sections.
%   Doesn't invoke FlowSim's own reader (that's exercised in unit tests).
fs_setup('smoke_mesh_load');

meshes = dir('*.msh');
fs_expect(numel(meshes) >= 10, ...
    sprintf('at least 10 .msh files at repo root (got %d)', numel(meshes)));

for k = 1:min(numel(meshes), 20)
    fn = meshes(k).name;
    try
        txt = fileread(fn);
        nNodes = numel(regexp(txt, '\$Nodes\s*\n(\d+)', 'once', 'tokens'));
        nElems = numel(regexp(txt, '\$Elements\s*\n(\d+)', 'once', 'tokens'));
        bytes = meshes(k).bytes;
        fprintf('  %-45s  %8d bytes\n', fn, bytes);
        fs_expect(bytes > 100, sprintf('%s > 100 bytes', fn));
    catch err
        fs_expect(false, sprintf('%s: read failed — %s', fn, err.message));
    end
end

fs_teardown();
