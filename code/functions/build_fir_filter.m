function kernel = build_fir_filter(sr, ged_params)
% Build FIR filter kernel for the GED target frequency band.

    nyquist = sr / 2;
    
    kernel = fir1(ged_params.filter_order, ...
        [ged_params.filter_lo / nyquist, ged_params.filter_hi / nyquist]);

end