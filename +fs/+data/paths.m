function p = paths(key)
%FS.DATA.PATHS  Resolve well-known FlowSim data file paths.
%
%   p = fs.data.paths('spe10')     % → path to spe10.mat if found, else ''
%   p = fs.data.paths('spe_perm')  % → path to spe_perm.dat if found, else ''
%
%   Search order:
%     1. FlowSim repo root (fsRoot appdata or pwd)
%     2. Owner's my-axon store (~/.../my-axon/dev-projects/flowsim-vectorize/
%        data/legacy-binaries/)
%     3. Environment variable FS_DATA_DIR
%
%   Owner-authorised Option D relocation (2026-07-03).
%   README at <my-axon>/dev-projects/flowsim-vectorize/data/legacy-binaries/README.md

    known = struct('spe10', 'spe10.mat', 'spe_perm', 'spe_perm.dat', 'gmsh', 'gmsh.exe');
    if ~isfield(known, key)
        error('fs.data.paths: unknown key ''%s'' (valid: %s)', ...
              key, strjoin(fieldnames(known), ', '));
    end
    filename = known.(key);

    candidates = {};

    fsRoot = getappdata(0, 'fs_test_root');
    if isempty(fsRoot), fsRoot = pwd; end
    candidates{end+1} = fullfile(fsRoot, filename); %#ok<AGROW>

    homeDir = getenv('HOME');
    if isempty(homeDir), homeDir = getenv('USERPROFILE'); end
    if ~isempty(homeDir)
        candidates{end+1} = fullfile(homeDir, 'projects', 'new-axon', 'axon', ...
            'my-axon', 'dev-projects', 'flowsim-vectorize', 'data', ...
            'legacy-binaries', filename); %#ok<AGROW>
    end

    envDir = getenv('FS_DATA_DIR');
    if ~isempty(envDir), candidates{end+1} = fullfile(envDir, filename); end %#ok<AGROW>

    for k = 1:numel(candidates)
        try
            if ~isempty(dir(candidates{k}))
                p = candidates{k};
                return;
            end
        catch
        end
    end

    p = '';
end
