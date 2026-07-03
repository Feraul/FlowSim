%UNIT_FS_ITER_WRAPPERS  Verify PR-F2 cross-cluster shim wrappers resolve.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_fs_iter_wrappers');

for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

wrappers = {'fs.iter.picard', 'fs.iter.anderson', 'fs.iter.lscheme', 'fs.lpew.dmpWeights'};
for k = 1:numel(wrappers)
    name = wrappers{k};
    w = which(name);
    fs_expect(~isempty(w), sprintf('%s exists on path', name));
end

fs_teardown();
