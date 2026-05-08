function peak_frequencies = retain_frequency_band(cho_output, retain_band)
% Retain CHO centre frequencies within a specified frequency band.
%
% Parameters
% ----------
% cho_output : struct
%   Output from run_cho_frequency_estimation.
%
% retain_band : [lo hi]
%   Frequency band to retain, in Hz.
%
% Returns
% -------
% peak_frequencies : cell array
%   Same size as cho_output.results.
%   Each cell contains valid centre frequencies in Hz, or NaN if no valid
%   oscillation was detected.
    
    cho_results = cho_output.results;
    
    peak_frequencies = cell(size(cho_results));
    
    for win_idx = 1:numel(cho_results)
    
        this_cho = cho_results{win_idx};
    
        % No CHO output or no detected oscillatory episodes
        if ~isfield(this_cho, 'bounding_boxes') || isempty(this_cho.bounding_boxes)
            peak_frequencies{win_idx} = NaN;
            continue
        end
    
        % Extract CHO centre frequencies
        all_freqs = [this_cho.bounding_boxes.center_fp]';
    
        % Retain frequencies within requested band
        valid_freqs = all_freqs( ...
            all_freqs >= retain_band(1) & ...
            all_freqs <= retain_band(2));
    
        if isempty(valid_freqs)
            peak_frequencies{win_idx} = NaN;
        else
            peak_frequencies{win_idx} = valid_freqs;
        end
    
    end

end