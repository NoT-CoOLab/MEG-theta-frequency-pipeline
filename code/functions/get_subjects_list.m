function subject_codes = get_subjects_list(paths_subfolders)
    subjDirs = dir(paths_subfolders);
    subjDirs = subjDirs([subjDirs.isdir]);
    
    % keep only subjects that have an MEG folder
    hasMEG = false(numel(subjDirs),1);
    for i = 1:numel(subjDirs)
        eegPath = fullfile(subjDirs(i).folder, subjDirs(i).name, 'meg');
        hasMEG(i) = exist(eegPath,'dir') == 7;
    end
    
    subject_codes = subjDirs(hasMEG); % e.g. ["101","102","103"...], or ["con1","con2","exp1"...]
    subject_codes = {subject_codes.name}';
    
    fprintf('Found %d subjects with MEG data in derivatives.\n', numel(subject_codes));