%UNIT_ASSEMBLY_SCAFFOLDS  Verify +fs/+assembly package scaffolds load.
%
%   PR-C4/C5/C6 scope: scaffold packages exist and are callable. Deep
%   parity tests (vs legacy + golden) will live in per-method unit tests
%   authored alongside each method's full vectorization PR (PR-C4b/5b/6b).

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_assembly_scaffolds');

for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

% Just verify the package functions resolve
scaffolds = {'fs.assembly.mpfah.build', ...
             'fs.assembly.nlfvpp.build', ...
             'fs.assembly.mpfaql.build', ...
             'fs.assembly.nlfvh.build', ...
             'fs.assembly.dmp.build'};
for k = 1:numel(scaffolds)
    name = scaffolds{k};
    try
        h = str2func(name);
        fs_expect(isa(h, 'function_handle'), sprintf('%s resolves as function_handle', name));
    catch err
        fs_expect(false, sprintf('%s failed to resolve: %s', name, err.message));
    end
end

fs_teardown();
