% run_theta_frequency_pipeline.m
%
% Hippocampal theta frequency estimation from MEG source-reconstructed data.
% Applies GED followed by CHO to estimate oscillatory centre frequencies
% within user-defined analysis windows.
%
% Before running:
%   1. Complete MEG source reconstruction (hippocampal_source_reconstruction.py)
%      to produce per-subject *_desc-lcmv-stc_*.mat files.
%   2. Add the functions/ folder to your MATLAB path.
%   3. Edit the defaults setup section.
%
% Output:
%   One .mat file per subject saved to output_dir, containing peak frequency
%   estimates per analysis window.
%
% 

clc; clear;
 
%% Setup paths

code_dir = 'path_to_the_functions'; % Contains this script and all functions called by it
addpath(genpath(code_dir)); % Add pipeline functions to path 
base_dir    = '/path_to_derivatives'; % Root derivatives folder
output_dir  = fullfile(base_dir, 'frequency_analysis'); % where to store the outputs
 
% Get the list of subject
subjects = get_subjects_list(base_dir);

% Ensure the main output folders exist
[~] = setup_dir(output_dir, false); % false = without a subfolder labelled as today's date

%% Setup filenames 

stc_pattern      = '%s_*desc-lcmv-stc*.mat';       % lcmv source time courses
metadata_pattern = '%s_*desc-lcmv-metadata*.mat';   % source space metadata
output_suffix    = 'desc-theta-freqs_meg.mat';       % output filename suffix

%% Setup parameters 

% Analysis window configuration
% GED and CHO are applied within analysis windows.
%
% Options:
%   'as_loaded' : use each loaded epoch directly as one analysis window
%   'fixed'     : subdivide each loaded epoch/continuous segment into
%                 fixed-length windows
%
% Examples:
%   Trial-based task:
%       cfg.window_mode = 'as_loaded';
%
%   Resting-state or pseudo real-time analysis:
%       cfg.window_mode = 'fixed';
%       cfg.window_size = 1000; % 1 second at 1000 Hz

cfg.window_mode = 'fixed';       % split each loaded epoch/continuous segment into fixed windows
cfg.window_size = 1000;          % samples, e.g. 1 second at 1000 Hz

% GED parameters
params.ged.filter_lo    = 2;    % high-pass edge of bandpass filter (Hz)
params.ged.filter_hi    = 8;    % low-pass edge of bandpass filter (Hz)
params.ged.filter_order = 1000; % FIR filter order 
params.ged.reg          = 0.01; % regularisation parameter for GED
params.ged.plot         = 0;    % 1 = plot covariance matrices and scree plot (disabled in chunked mode)

% CHO parameters
params.cho.minimum_cycles   = 2;
params.cho.ovlp_threshold   = 0.5;
params.cho.frequency_vector = 1:40;
params.cho.plot             = 0;

% Frequency band to retain from CHO output
cfg.retain_band = [2 8]; % Only peaks within this range are kept (Hz)
 
%% Run the hippocampal frequency estimation pipeline

for sub_idx = 1:length(subjects)
    subj = subjects{sub_idx};
    fprintf('\n[%d/%d] Processing %s ...\n', sub_idx, length(subjects), subj);
 
    % Locate input data
    input_folder = fullfile(base_dir, 'source_reconstruction', subj);
    stc_files  = dir(fullfile(input_folder, sprintf(stc_pattern, subj)));
    
    if isempty(stc_files)
        warning('[%s] No lcmv-stc file found in %s.. skipping.', subj, input_folder);
        continue
    end
 
    % Locate metadata file
    meta_files = dir(fullfile(input_folder, sprintf(metadata_pattern, subj)));
  
    if isempty(meta_files)
        warning('[%s] No metadata file found.. skipping.', subj);
        continue
    end
 
    % Subject-level config
    cfg.stc_path      = fullfile(input_folder, stc_files(1).name);
    cfg.metadata_path = fullfile(input_folder, meta_files(1).name);
    cfg.output_folder = fullfile(output_dir, subj);
    [~] = setup_dir(cfg.output_folder, false);

    % Run frequency estimation
    peak_frequencies = estimate_hippocampal_theta_frequency(cfg, params);
 
    % Save results
    outname = sprintf('%s_%s', subj, output_suffix);
    save(fullfile(cfg.output_folder, outname), 'peak_frequencies', 'cfg', 'params');

    fprintf('[%s] Done. Results saved to %s\n', subj, cfg.output_folder);
end
 
