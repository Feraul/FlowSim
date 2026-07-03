%UNIT_BASELINE_REPRODUCES  Rerun capture_baseline and diff against committed golden .mat.
%
%   The correctness oracle in one file: recapture the numerical fingerprint,
%   load the committed baseline, assert bit-identical structural fields +
%   Frobenius-tight numeric fields.
%
%   Fails if the current run diverges from the baseline. First run captures
%   the golden; subsequent runs verify reproducibility.

addpath(fullfile(pwd, 'tests', 'helpers'));
fs_setup('unit_baseline_reproduces');

setappdata(0, 'fs_test_root', pwd);
fsRoot = pwd;   % snapshot — fs_test_env will cd away, but goldens live under fsRoot
for sub = {'','base','solvers','factories','simulacoes','benchmarks'}
    d = fullfile(pwd, sub{1});
    if isfolder(d), addpath(d); end
end

meshes = {'M8.msh'};
methods = {'tpfa', 'mpfad'};
numcase = 439;

for mi = 1:numel(meshes)
    for mj = 1:numel(methods)
        mesh = meshes{mi}; method = methods{mj};
        stem = sprintf('%s-num%g-%s', strrep(mesh, '.msh', ''), numcase, method);
        goldenFile = fullfile(fsRoot, 'tests', 'golden', [stem '.mat']);

        % Try to load — bypasses WSL/UNC exist()/dir()/isfile() cache issues.
        % A successful load proves the file is readable; a caught error means missing.
        try
            gold = load(goldenFile);
            gold = gold.capture;
        catch err
            fs_expect(false, sprintf('golden not loadable: %s (%s)', stem, err.message));
            continue;
        end

        % Recapture into a scratch outfile (don't overwrite golden)
        tmpFile = fullfile(tempdir, [stem '.tmp.mat']);
        [~, cap] = capture_baseline('mesh', mesh, 'method', method, ...
                                    'numcase', numcase, 'outfile', tmpFile);

        % Structural equality
        fs_expect(cap.nNodes  == gold.nNodes,  sprintf('%s: nNodes matches (%d)',  stem, cap.nNodes));
        fs_expect(cap.nElems  == gold.nElems,  sprintf('%s: nElems matches (%d)',  stem, cap.nElems));
        fs_expect(cap.nBFaces == gold.nBFaces, sprintf('%s: nBFaces matches (%d)', stem, cap.nBFaces));
        fs_expect(cap.nIFaces == gold.nIFaces, sprintf('%s: nIFaces matches (%d)', stem, cap.nIFaces));
        fs_expect(cap.nCorners == gold.nCorners, sprintf('%s: nCorners matches (%d)', stem, cap.nCorners));

        % Assembly fingerprint (Frobenius diff)
        if isfield(gold, 'M_frobnorm') && isfield(cap, 'M_frobnorm')
            fs_expect(cap.M_nnz == gold.M_nnz, ...
                sprintf('%s: M nnz matches (%d)', stem, cap.M_nnz));

            frobRel = abs(cap.M_frobnorm - gold.M_frobnorm) / max(1e-30, gold.M_frobnorm);
            fs_expect(frobRel < 1e-12, ...
                sprintf('%s: M Frobenius reproducible (rel diff %.3e)', stem, frobRel));

            iRel = abs(cap.I_norm - gold.I_norm) / max(1e-30, gold.I_norm);
            fs_expect(iRel < 1e-10, ...
                sprintf('%s: I L2-norm reproducible (rel diff %.3e)', stem, iRel));
        end

        % Premethod per-field norms
        if isfield(gold, 'premethod_frobs') && isfield(cap, 'premethod_frobs')
            fs_expect(numel(gold.premethod_frobs) == numel(cap.premethod_frobs), ...
                sprintf('%s: premethod field count matches', stem));
            for k = 1:min(numel(gold.premethod_frobs), numel(cap.premethod_frobs))
                gv = gold.premethod_frobs{k};
                cv = cap.premethod_frobs{k};
                if isnan(gv) && isnan(cv), continue; end
                if isnan(gv) || isnan(cv)
                    fs_expect(false, sprintf('%s: premethod.%s nan-mismatch', stem, gold.premethod_keys{k}));
                    continue;
                end
                rel = abs(gv - cv) / max(1e-30, gv);
                fs_expect(rel < 1e-10, ...
                    sprintf('%s: premethod.%s L2 reproducible (rel %.3e)', stem, gold.premethod_keys{k}, rel));
            end
        end

        % Cleanup tmp
        if isfile(tmpFile), delete(tmpFile); end
    end
end

fs_teardown();
