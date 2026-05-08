function ged_output = run_ged(analysis_windows, ged_params)
% apply_ged
%
% Applies generalized eigendecomposition (GED) to identify the spatial
% component that maximally expresses theta-band activity.
%
% The GED contrasts:
%   S = covariance of theta-band filtered hippocampal source activity
%   R = covariance of broadband hippocampal source activity
%
% This implementation follows the general GED framework used for
% oscillatory component enhancement in electrophysiological data, as
% described in Mike X Cohen's teaching materials and related published work.
% It is adapted here for hippocampal source-reconstructed MEG data arranged
% as sources x time within each analysis window.
%
% Inputs
% ------
% filtered_signal  : hippocampal sources x time
%   Theta-band filtered source activity used to construct covariance S.
%
% broadband_signal : hippocampal sources x time
%   Broadband source activity used to construct covariance R.
%
% reg : scalar
%   Regularisation parameter applied to R.
%
% do_plot : logical
%   If true, plot covariance matrices and GED eigenvalue spectrum.
%
% Output
% ------
% theta_component : 1 x time
%   GED-enhanced theta component time course.

sr = analysis_windows.sampling_rate;

% Build theta-band FIR filter kernel
kernel = build_fir_filter(sr, ged_params);

window_data = analysis_windows.data;
ged_components = cell(size(window_data));

% Avoid too many figures in fixed-window mode
if strcmp(analysis_windows.mode, 'fixed') && ged_params.plot
    warning('[GED] Plotting is disabled in fixed-window mode. It would create too many figures.');
    ged_params.plot = 0;
end

for win_idx = 1:numel(window_data)

    signal = window_data{win_idx};

    % Filtered data define the theta-band covariance
    filtered_signal = apply_fir_filter(signal, kernel);

    % GED extracts the component that maximally expresses theta activity
    ged_components{win_idx} = apply_ged( ...
        filtered_signal, ...
        signal, ...
        ged_params.reg, ...
        ged_params.plot);

end

ged_output.components = ged_components;
ged_output.sampling_rate = sr;
ged_output.mode = analysis_windows.mode;
ged_output.window_size = analysis_windows.window_size;

end