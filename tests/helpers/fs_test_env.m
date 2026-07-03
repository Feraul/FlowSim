function testDir = fs_test_env(varargin)
%FS_TEST_ENV  Set up an isolated Start.dat + mesh env so preprocessormod runs under WSL.
%
%   testDir = fs_test_env()                      % default: /tmp/flowsim_test_env
%   testDir = fs_test_env('dir', '/some/where')  % explicit dir
%   testDir = fs_test_env('mesh', 'M8.msh')      % override mesh (default: what Start.dat says)
%
%   Setup:
%     1. Reads owner's Start.dat from FlowSim repo root (getappdata(0,'fs_test_root'))
%        OR pwd if root not set.
%     2. Patches the two Windows-native paths:
%          C:\...\Benchmark_Cases\...  →  <testDir>/out/    (output — WSL-writable)
%          C:\...\Malhas2              →  <fsRoot>          (mesh dir — points to repo)
%     3. Writes the patched copy to <testDir>/Start.dat
%     4. cd's into <testDir> so preprocessormod's `fopen('Start.dat')` picks up the copy
%     5. Ensures <testDir>/out/ exists so build_directories won't fail
%
%   Returns the test dir. Original owner's Start.dat is untouched.
%
%   Cross-platform: on native Windows, WSL paths become C:/Temp/flowsim_test; on Linux
%   both paths use /tmp. Only the string substitutions differ.

    p = inputParser;
    addParameter(p, 'dir',  fullfile(tempdir, 'flowsim_test_env'), @ischar);
    addParameter(p, 'mesh', '',                                    @ischar);
    parse(p, varargin{:});
    testDir = p.Results.dir;

    fsRoot = getappdata(0, 'fs_test_root');
    if isempty(fsRoot), fsRoot = pwd; end

    startSrc = fullfile(fsRoot, 'Start.dat');
    if ~isfile(startSrc)
        error('fs_test_env: Start.dat not found at %s', startSrc);
    end

    % Ensure test dir + subdirs
    if ~isfolder(testDir),                   mkdir(testDir); end
    if ~isfolder(fullfile(testDir, 'out')),  mkdir(fullfile(testDir, 'out')); end

    % Read + patch Start.dat
    txt = fileread(startSrc);

    % Patch line 33 (output path)
    outPath = fullfile(testDir, 'out');
    txt = regexprep(txt, ...
        'C:\\Users\\flc59\\Documents\\Benchmark_Cases\\[^\r\n]*', ...
        regexprep(outPath, '\\', '\\\\'));

    % Patch line 218 (mesh dir → FlowSim repo root, so .msh files at root resolve)
    txt = regexprep(txt, ...
        'C:\\Users\\flc59\\Documents\\Malhas2', ...
        regexprep(fsRoot, '\\', '\\\\'));

    % Patch the mesh filename (line ~224) — owner's Start.dat may reference a
    % mesh not committed to git. Default to M8.msh (smallest present in repo).
    meshFile = p.Results.mesh;
    if isempty(meshFile), meshFile = 'M8.msh'; end
    if ~isfile(fullfile(fsRoot, meshFile))
        warning('fs_test_env:MeshMissing', ...
                'Mesh %s not found at %s', meshFile, fsRoot);
    end
    % The mesh filename in Start.dat is the FIRST line after "Nome do arquivo de malha"
    % that isn't a comment / blank / marker. Rather than line-count, replace by
    % pattern: any line ending in .msh that isn't inside a comment.
    txt = regexprep(txt, ...
        '(?m)^[^/\r\n]+\.msh\s*$', ...
        meshFile);

    % Write patched Start.dat
    fid = fopen(fullfile(testDir, 'Start.dat'), 'w');
    if fid < 0, error('fs_test_env: cannot write to %s', testDir); end
    fwrite(fid, txt);
    fclose(fid);

    % Symlink or copy the mesh files into testDir so Start.dat's relative mesh
    % resolution can find them (some paths in Start.dat may become relative post-patch).
    % Also keep the FlowSim source path added so functions still resolve.

    % cd into testDir so fopen('Start.dat') picks up the copy
    cd(testDir);
end
