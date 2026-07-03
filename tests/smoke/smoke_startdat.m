addpath(fullfile(pwd, 'tests', 'helpers'));  %% bootstrap so fs_* helpers resolve
%SMOKE_STARTDAT  Verify Start.dat parses and report key config values.
%
%   Uses whatever parser FlowSim's own preprocessormod uses (via getdatafile).
%   This is the fix for study gap #1 (unknown pmethod).
fs_setup('smoke_startdat');

fs_expect(exist('Start.dat', 'file') == 2, 'Start.dat exists at repo root');
fs_expect(exist('getdatafile', 'file') == 2, 'getdatafile function on path');

% Try to invoke getdatafile — signature varies; try common calls
try
    cfg = getdatafile(0);
    fs_expect(true, 'getdatafile(0) returned');

    if isstruct(cfg)
        if isfield(cfg, 'numcase'), fprintf('  numcase   = %g\n', cfg.numcase); end
        if isfield(cfg, 'pmethod'), fprintf('  pmethod   = %s\n', char(cfg.pmethod)); end
        if isfield(cfg, 'phasekey'), fprintf('  phasekey  = %g\n', cfg.phasekey); end
        if isfield(cfg, 'smethod'), fprintf('  smethod   = %s\n', char(cfg.smethod)); end
        if isfield(cfg, 'acel'), fprintf('  acel      = %s\n', char(cfg.acel)); end
    else
        fprintf('  getdatafile returned non-struct: %s\n', class(cfg));
    end
catch err
    fs_expect(false, sprintf('getdatafile(0) failed — %s', err.message));
end

fs_teardown();
