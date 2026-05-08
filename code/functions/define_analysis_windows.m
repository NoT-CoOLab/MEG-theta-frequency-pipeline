function analysis_windows = define_analysis_windows(source_data, cfg)
% Define the analysis windows used for GED and CHO.
%
% Window modes
% ------------
% 'as_loaded'
%   Uses each loaded epoch directly as one analysis window.
%
% 'fixed'
%   Splits each loaded epoch (or continuous segment) into fixed-length
%   windows of cfg.window_size samples.
%
% Returns
% -------
% analysis_windows : struct
%   .data           : cell array containing analysis windows
%   .mode           : analysis window mode
%   .window_size    : window size in samples (fixed mode only)
%   .sampling_rate  : sampling rate in Hz

epochs = source_data.epochs;
n_epochs = numel(epochs);

switch lower(cfg.window_mode)

    % Use loaded epochs directly
    case 'as_loaded'

        window_data = cell(n_epochs, 1);

        for epoch_idx = 1:n_epochs
            window_data{epoch_idx, 1} = epochs{epoch_idx};
        end

        analysis_windows.mode = 'as_loaded';
        analysis_windows.window_size = [];

    
    % Split data into fixed-length windows
    case 'fixed'

        window_size = cfg.window_size;

        n_windows = floor(size(epochs{1}, 2) / window_size);

        fprintf('[Windows] Fixed-window mode: %d windows of %d samples per epoch\n', ...
            n_windows, window_size);

        window_data = cell(n_epochs, n_windows);

        for epoch_idx = 1:n_epochs

            for window_idx = 1:n_windows

                t_start = (window_idx - 1) * window_size + 1;
                t_end   =  window_idx      * window_size;

                window_data{epoch_idx, window_idx} = ...
                    epochs{epoch_idx}(:, t_start:t_end);

            end
        end

        analysis_windows.mode = 'fixed';
        analysis_windows.window_size = window_size;

    otherwise

        error(['Unknown cfg.window_mode: %s\n' ...
               'Valid options are ''as_loaded'' or ''fixed''.'], ...
               cfg.window_mode);

end

analysis_windows.data = window_data;
analysis_windows.sampling_rate = source_data.sampling_rate;

end