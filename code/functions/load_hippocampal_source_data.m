function source_data = load_hippocampal_source_data(cfg)

    loaded_stc  = load(cfg.stc_path);
    loaded_meta = load(cfg.metadata_path);
    metadata    = loaded_meta.metadata;
    
    epoch_fields = fieldnames(loaded_stc.epochs);
    n_epochs = numel(epoch_fields);
    
    epochs = cell(n_epochs, 1);
    
    for epoch_idx = 1:n_epochs
        epochs{epoch_idx} = loaded_stc.epochs.(epoch_fields{epoch_idx});
    end
    
    source_data.epochs = epochs;
    source_data.metadata = metadata;
    source_data.sampling_rate = metadata.sampling_rate;
    source_data.chunk_size = cfg.chunk_size;

end