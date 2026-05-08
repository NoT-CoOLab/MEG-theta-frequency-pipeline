function [outputs] = oof_fitting_v2(psd_in_db,frequency_domain)
% FOOOF_MATLAB_V1 
% oof_fitting_v2: This function fits a power spectral density (PSD) signal to a model 
% in the frequency domain by separating aperiodic components from periodic components.
% psd_in_db and frequency_domain should be column vector
% frequency_domain should avoid 0

%% functions for fitting
modelfun = @(b,x) b(1) - log10(b(2) + x(:,1).^b(3));

% This defines the model function, which is a log-linear function of frequency.
% The model is parameterized by three variables: b(1), b(2), and b(3).
% b(1): Amplitude offset (intercept in log space)
% b(2): Baseline offset that prevents log of zero
% b(3): Exponent determining the slope

%% settings
gaussians = [];
outputs = [];
options = optimoptions('lsqcurvefit','Display','off');
MAX_GAUSSIANS = 4; %unused
X = frequency_domain;
Y = psd_in_db;
ap_percentile_thresh = 10; % percentile threshold for aperiodic component separation

%% Fit aperiodic signal - crucial and not trivial
beta0 = [min(Y) 0 1]; % [Amplitude offset, Baseline, Slope exponent]
lb = [min(Y)-2 0 0]; % Lower bounds for the parameters
ub = [max(Y)+2 0 1000 ]; % Upper bounds for the parameters

% perform initial curve fitting using nonlinear least squares
mdl_coef_init = lsqcurvefit(modelfun,beta0,X,Y,lb,ub,options);
initial_fit = modelfun(mdl_coef_init, X);

% compute the "flattened" spectrum by subtracting the aperiodic fit from the original PSD
flatspec = Y - initial_fit;

% remove negative values from the flattened spectrum (as they don't make sense for power)
flatspec(flatspec<0) = 0;
% identify a threshold for ignoring aperiodic components based on a percentile
perc_thresh = prctile(flatspec,ap_percentile_thresh);
% mask to ignore certain frequencies where the flattened spectrum is below the threshold
perc_mask = flatspec <= perc_thresh;
freqs_ignore = X(perc_mask); % Frequencies to ignore
spectrum_ignore = Y(perc_mask); % Corresponding PSD values to ignore

% If there are enough points, refit the model to the filtered data
if length(spectrum_ignore) > 3

    mdl_coef = lsqcurvefit(modelfun,beta0,freqs_ignore,spectrum_ignore,lb,ub,options);
    initial_fit2 = modelfun(mdl_coef, X);
    
    %% results
    outputs.initial_ap_fit_coef = mdl_coef; % Store the coefficients of the fit
    outputs.init_ap_fit = initial_fit2; % Store the fit result
    % outputs.final_fit = final_fit;
else
    outputs.initial_ap_fit_coef = nan;% If not enough data points, return NaN values
    outputs.init_ap_fit = NaN(size(psd_in_db));
end

%% visualise result
% figure();
% plot(X,psd_in_db);
% hold on;
% plot(X,initial_fit2);
% fprintf('error rate = %d\n',error_rate);

end

