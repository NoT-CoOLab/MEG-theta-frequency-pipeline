function cho_output = run_cho(ged_output, cho_params)
% Apply CHO to each GED-enhanced theta component.
%
% CHO is used to detect genuine oscillations within each analysis window and
% estimate their centre frequency.

    components = ged_output.components;
    sr = ged_output.sampling_rate;
    
    cho_results = cell(size(components));
    
    % Avoid too many figures in fixed-window mode
    if strcmp(ged_output.mode, 'fixed') && cho_params.plot
        warning('[CHO] Plotting is disabled in fixed-window mode.');
        cho_params.plot = 0;
    end
    
    for win_idx = 1:numel(components)
    
        theta_component = components{win_idx};
    
        try
            cho_results{win_idx} = CHO_v22_nr(theta_component, sr, cho_params);
        catch ME
            warning('[CHO] CHO failed for window %d: %s', win_idx, ME.message);
            cho_results{win_idx}.bounding_boxes = [];
        end
    
    end
    
    cho_output.results = cho_results;
    cho_output.sampling_rate = sr;
    cho_output.mode = ged_output.mode;
    cho_output.window_size = ged_output.window_size;

end