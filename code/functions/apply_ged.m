function theta_component = apply_ged(filtered_signal, broadband_signal, reg, do_plot)
% apply_ged
%
% Applies generalized eigendecomposition (GED) to identify the spatial
% component that maximally expresses theta-band activity.
%
% GED contrasts:
%   S = covariance of theta-band filtered activity
%   R = covariance of broadband activity
%
% Inputs
% ------
% filtered_signal : matrix (sources x time)
%   Theta-band filtered hippocampal source activity.
%
% broadband_signal : matrix (sources x time)
%   Broadband hippocampal source activity.
%
% reg : scalar
%   Regularisation parameter for covariance matrix R.
%
% do_plot : logical
%   If true, plot covariance matrices and GED eigenspectrum.
%
% Output
% ------
% theta_component : vector (time x 1)
%   GED-enhanced theta component time course.

    %% Defaults
    
    if nargin < 4
        do_plot = false;
    end
    
    if nargin < 3 || isempty(reg)
        reg = 0.01;
    end
    
    %% Construct covariance matrix S (theta-band filtered activity)
    
    covTheta = cov(filtered_signal');
    
    %% Construct covariance matrix R (broadband activity)
    
    covBroadband = cov(broadband_signal');
    
    %% Regularize covariance matrix R
    
    covR = (1 - reg) * covBroadband + ...
            reg * mean(eig(covBroadband)) * eye(size(filtered_signal,1));
    
    %% Perform generalized eigendecomposition
    
    [eigenvecs, eigenvalues] = eig(covTheta, covR);
    
    % Sort eigenvalues in descending order
    [eigenvalues_sorted, sort_idx] = sort(diag(eigenvalues), 'descend');
    
    eigenvecs = eigenvecs(:, sort_idx);
    
    %% Normalize eigenvectors to unit length
    
    eigenvecs = bsxfun(@rdivide, ...
                       eigenvecs, ...
                       sqrt(sum(eigenvecs.^2, 1)));
    
    %% Optional plotting
    
    if do_plot
    
        figure('Name', 'GED diagnostics');
    
        subplot(1,3,1)
        imagesc(covTheta)
        axis square
        xlabel('Source')
        ylabel('Source')
        title('Theta covariance (S)')
        colorbar
    
        subplot(1,3,2)
        imagesc(covR)
        axis square
        xlabel('Source')
        ylabel('Source')
        title('Broadband covariance (R)')
        colorbar
    
        subplot(1,3,3)
        plot(eigenvalues_sorted, 'ko-', ...
             'MarkerFaceColor', 'w', ...
             'LineWidth', 1.5)
    
        xlabel('Component')
        ylabel('\lambda')
        title('GED eigenspectrum')
    
    end
    
    %% Compute dominant theta component time series
    
    theta_component = eigenvecs(:,1)' * broadband_signal;
    
    % Return as column vector (time x 1)
    theta_component = theta_component';

end