function peak_frequencies = estimate_hippocampal_theta_frequency(cfg, params)

    % 1. Load hippocampal source data
    source_data = load_hippocampal_source_data(cfg);

    % 2. Split data into analysis segments: full epochs or chunks
    analysis_windows = define_analysis_windows(source_data, cfg);
    
    % 3. Apply GED to enhance theta activity in each segment
    ged_output = run_ged(analysis_windows, params.ged);
    
    % 4. Apply CHO to detect genuine oscillations and estimate centre frequency
    cho_output = run_cho(ged_output, params.cho);
    
    % 5. Retain valid peaks in frequency band of interest
    peak_frequencies = retain_frequency_band(cho_output, cfg.retain_band);

end