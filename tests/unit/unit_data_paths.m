%UNIT_DATA_PATHS  Verify fs.data.paths resolves known keys and rejects unknown.
addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_data_paths');

% Known keys should not error
for k = {'spe10', 'spe_perm', 'gmsh'}
    try
        p = fs.data.paths(k{1});
        fs_expect(true, sprintf('fs.data.paths(''%s'') returns without error', k{1}));
        if ~isempty(p)
            fprintf('  (resolved) %s → %s\n', k{1}, p);
        else
            fprintf('  (not found) %s\n', k{1});
        end
    catch err
        fs_expect(false, sprintf('%s errored: %s', k{1}, err.message));
    end
end

% Unknown key should error clearly
gotErr = false;
try
    fs.data.paths('nonexistent_key');
catch
    gotErr = true;
end
fs_expect(gotErr, 'unknown key triggers assertion');

fs_teardown();
