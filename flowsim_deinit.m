function flowsim_deinit(varargin)
%FLOWSIM_DEINIT  Remove all FlowSim source folders from the MATLAB path.
%
%   Symmetric to flowsim_init. Reads the paths stashed in appdata by init
%   and rmpath's each. Safe to call multiple times (no-op if not init'd).
%
%   flowsim_deinit
%   flowsim_deinit('verbose', true)

    p = inputParser;
    addParameter(p, 'verbose', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    paths = getappdata(0, 'flowsim_init_paths');
    if isempty(paths)
        if p.Results.verbose
            fprintf('flowsim_deinit: no init state found (no-op).\n');
        end
        return;
    end

    for k = 1:numel(paths)
        try
            rmpath(paths{k});
            if p.Results.verbose
                fprintf('  - %s\n', paths{k});
            end
        catch err
            warning('flowsim_deinit:RmpathFailed', ...
                    'Could not remove %s: %s', paths{k}, err.message);
        end
    end

    rmappdata(0, 'flowsim_init_paths');
    rmappdata(0, 'flowsim_init_root');

    if p.Results.verbose
        fprintf('flowsim_deinit: cleaned %d path entries.\n', numel(paths));
    end
end
