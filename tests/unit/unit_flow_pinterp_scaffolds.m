%UNIT_FLOW_PINTERP_SCAFFOLDS  Verify Phase D scaffolds resolve.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_flow_pinterp_scaffolds');

for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

% Existence check via which — more reliable than str2func
fns = {'fs.flow.mpfad', 'fs.flow.tpfa', 'fs.lpew.pinterp'};
for k = 1:numel(fns)
    name = fns{k};
    w = which(name);
    fs_expect(~isempty(w), sprintf('%s exists on path (%s)', name, w));
end

fs_teardown();
