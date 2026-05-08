# -*- coding: utf-8 -*-
"""
hippocampal_source_reconstruction.py

Performs MEG source reconstruction targeting hippocampal sources using a
mixed source space (bilateral hippocampi + sparse cortical surface) and
an LCMV beamformer. Outputs epoch-wise hippocampal source time series as
.mat files for downstream analysis in MATLAB (GED + CHO).

Pipeline steps (each step checks for existing output and skips if found):
    1. BEM solution
    2. Coregistration (launches GUI if trans file not found)
    3. Mixed source space (hippocampi + sparse cortex)
    4. Forward model
    5. LCMV beamformer + save hippocampal source time series as .mat

Requirements:
    - FreeSurfer reconstruction completed for each subject
    - Preprocessed MEG epochs (.fif) in BIDS-compatible folder structure
    - MNE-Python >= 1.6 (see environment.yml)

Usage:
    Configure the CONFIGURATION section below, then run:
        python hippocampal_source_reconstruction.py

@author: eleonorm
"""

#%% CONFIGURATION
# Root directories
MRI_DIR             = r'/path/to/mri_data'          # FreeSurfer subjects dir
PREPROC_DATA_PATH   = r'/path/to/derivatives/preprocessing'
SOURCE_RESULTS_PATH = r'/path/to/derivatives/source_reconstruction'
RAW_DATA_PATH       = r'/path/to/raw_data_bids'

# Empty room recording path for noise covariance estimation.
# Set to the full path of the empty room .fif file for this subject/session.
# If None, an ad hoc noise covariance will be used instead.
EMPTY_ROOM_PATH = None

# Subject / session identifiers
SUBJ = '01'           # subject ID
TASK = 'rest'   # task label in BIDS filename 
RUN  = '001'           # run number (set to None if your files have no run label)

# Preprocessed epochs filename descriptor
EPOCHS_DESC = 'denoised_meg'

# Hippocampal labels to include in the volumetric source space
HIPPOCAMPAL_LABELS = ['Left-Hippocampus', 'Right-Hippocampus']

# Cortical surface source space density ('oct3' = sparse)
SURFACE_SPACING = 'oct3'

# BEM conductivity (single layer is sufficient for MEG)
BEM_CONDUCTIVITY = (0.3,)

# LCMV parameters
LCMV_REG        = 0.05   # regularisation parameter
LCMV_N_JOBS     = 4      # parallel jobs for covariance computation

#%% Imports
import os
import re
import numpy as np
import mne
from mne import setup_source_space, setup_volume_source_space
from mne.cov import compute_covariance
from mne.beamformer import make_lcmv, apply_lcmv_epochs
from scipy.io import savemat

#%% Helpers
def get_bids_filename(subj, task, run, desc, suffix='epo', ext='.fif'):
    """Build a BIDS-compatible filename from components."""
    run_str = f"_run-{run}" if run else ""
    return f"sub-{subj}_task-{task}{run_str}_desc-{desc}_{suffix}{ext}"

def get_source_dir(source_results_path, subj):
    """Return the subject-level source reconstruction directory, creating it if needed."""
    source_dir = os.path.join(source_results_path, f"sub-{subj}")
    os.makedirs(source_dir, exist_ok=True)
    return source_dir

#%% Create the forward model
"""
To compute the forward operator we need:
    1. the BEM surfaces
    2. a _trans.fif file that contains the coregistration info
    3. a source space
"""
# compute BEM solution
def compute_or_load_bem(subj, mri_dir, source_dir):
    """
    Compute the BEM solution from FreeSurfer surfaces, or load it if already saved.
    A single-layer BEM (inner skull only) is used, which is standard for MEG.
    Saves a plot of the BEM surfaces.
    """
    bids_subj  = f"sub-{subj}"
    mri_subj   = os.path.join(mri_dir, subj + '_converted')
    bem_path   = os.path.join(source_dir, f"{bids_subj}_bem.fif")

    if os.path.exists(bem_path):
        print(f"[BEM] Loading existing solution from {bem_path}")
        bem = mne.read_bem_solution(bem_path)
    else:
        print("[BEM] Computing BEM solution ...")
        model = mne.make_bem_model(
            mri_subj, ico=4,
            conductivity=BEM_CONDUCTIVITY,
            subjects_dir=mri_dir
        )
        bem = mne.make_bem_solution(model)
        mne.write_bem_solution(bem_path, bem, overwrite=True)

        # plot
        fig = mne.viz.plot_bem(
            subject=mri_subj, subjects_dir=mri_dir,
            brain_surfaces='white', orientation='coronal'
        )
        fig.savefig(
            os.path.join(source_dir, f"{bids_subj}_bem_coronal.png"),
            dpi=300, bbox_inches='tight'
        )
        print(f"[BEM] Solution saved to {bem_path}")

    return bem, mri_subj 

# Co-registration
def load_or_run_coregistration(subj, mri_subj, mri_dir, source_dir, raw_fif_path):
    """
    Load the coregistration transform (.fif), or launch the MNE coregistration
    GUI if it has not been run yet. The GUI must be completed manually before
    the script can continue.

    Note: coregistration is a manual step. The trans file should be saved to
    source_dir before re-running this script.
    """
    bids_subj = f"sub-{subj}"
    trans_path = os.path.join(source_dir, f"{bids_subj}_trans.fif")

    if not os.path.exists(trans_path):
        print("[COREG] No trans file found. Launching coregistration GUI ...")
        print(f"        Save the trans file to: {trans_path}")
        mne.utils.set_config("SUBJECTS_DIR", mri_dir, set_env=True)
        mne.gui.coregistration()
        
    else:
        # allow for any trans file in folder (e.g. with different naming)
        trans_files = [f for f in os.listdir(source_dir) if f.endswith('trans.fif')]
        trans_path = os.path.join(source_dir, trans_files[0])
        print(f"[COREG] Trans file loaded from {trans_path}")

    # plot (verify sensor/MRI alignment)
    info = mne.io.read_info(raw_fif_path)
    fig = mne.viz.plot_alignment(
        info, trans_path, subject=mri_subj, dig=True,
        meg=['helmet', 'sensors'], subjects_dir=mri_dir,
        surfaces='head-dense'
    )
    fig.plotter.screenshot(os.path.join(source_dir, f"{bids_subj}_coreg.png"))
    fig.plotter.close()
    print(f"[COREG] Alignment QC saved to {source_dir}")

    return trans_path

# Create a mixed source space
def compute_or_load_source_space(subj, mri_subj, mri_dir, source_dir, bem):
    """
    Build a mixed source space combining:
      - Bilateral hippocampal volumetric sources (from FreeSurfer aseg)
      - Sparse cortical surface sources (oct3)
    Or load from file if already computed.
    """
    bids_subj    = f"sub-{subj}"
    src_path     = os.path.join(source_dir, f"{bids_subj}_src.fif")
 
    if os.path.exists(src_path):
        print(f"[SOURCE SPACE] Loading existing source space from {src_path}")
        mixed_src = mne.read_source_spaces(src_path)
    else:
        print("[SOURCE SPACE] Computing mixed source space ...")
 
        # Volumetric hippocampal sources
        hip_src = setup_volume_source_space(
            subject=mri_subj, mri='aseg.mgz', bem=bem,
            volume_label=HIPPOCAMPAL_LABELS,
            subjects_dir=mri_dir
        )
 
        # Sparse cortical surface sources
        surf_src = setup_source_space(
            subject=mri_subj, subjects_dir=mri_dir,
            spacing=SURFACE_SPACING, add_dist=False
        )
 
        # Combine: cortical first, hippocampal appended
        mixed_src  = surf_src
        mixed_src += hip_src
 
        # plot
        fig = mne.viz.plot_alignment(
            subject=mri_subj, subjects_dir=mri_dir,
            surfaces=dict(white=0.3), coord_frame='mri',
            src=mixed_src
        )
        mne.viz.set_3d_view(fig, azimuth=180, elevation=30, distance=0.40)
        fig.plotter.screenshot(
            os.path.join(source_dir, f"{bids_subj}_src.png")
        )
 
        mne.write_source_spaces(src_path, mixed_src, overwrite=False)
        print(f"[SOURCE SPACE] Saved to {src_path}")
 
    return mixed_src

# Create the forward model
def compute_or_load_forward(subj, mri_subj, source_dir, raw_fif_path,
                             trans_path, mixed_src, bem):
    """
    Compute the forward solution (lead-field matrix) for the mixed source space,
    or load it if already computed.
    """
    bids_subj = f"sub-{subj}"
    fwd_path  = os.path.join(source_dir, f"{bids_subj}_fwd.fif")
 
    if os.path.exists(fwd_path):
        print(f"[FORWARD] Loading existing forward model from {fwd_path}")
        fwd = mne.read_forward_solution(fwd_path)
    else:
        print("[FORWARD] Computing forward model ...")
        fwd = mne.make_forward_solution(
            raw_fif_path,
            trans=trans_path,
            src=mixed_src,
            bem=bem,
            meg=True, eeg=False,
            mindist=5.,    # minimum distance (mm) from inner skull surface
            n_jobs=LCMV_N_JOBS,
            verbose=True
        )
        mne.write_forward_solution(fwd_path, fwd, overwrite=False)
        print(f"[FORWARD] Saved to {fwd_path}")
 
    return fwd

#%% LCMV beamformer

# noise covariance
def compute_noise_cov(epochs_meg, empty_room_path):
    """
    Compute noise covariance from an empty room recording if provided,
    falling back to an ad hoc covariance otherwise.
 
    Parameters
    ----------
    epochs_meg : mne.Epochs
        Preprocessed MEG epochs (used for channel info and ad hoc fallback).
    empty_room_path : str or None
        Full path to the empty room .fif file for this subject/session.
        Set EMPTY_ROOM_PATH in the user configuration block. If None,
        ad hoc covariance is used.
    """
    if empty_room_path is not None and os.path.exists(empty_room_path):
        print(f"[NOISE COV] Loading empty room recording from {empty_room_path}")
        empty_room_raw = mne.io.read_raw_fif(empty_room_path, preload=True)
        # Match channel selection to MEG epochs
        empty_room_raw.pick(epochs_meg.ch_names)
        noise_cov = mne.compute_raw_covariance(
            empty_room_raw, method='shrunk', rank='info'
        )
        print("[NOISE COV] Noise covariance computed from empty room recording.")
    else:
        if empty_room_path is not None:
            print(f"[NOISE COV] Warning: empty room file not found at {empty_room_path}")
        print("[NOISE COV] Using ad hoc noise covariance.")
        noise_cov = mne.make_ad_hoc_cov(epochs_meg.info)
 
    return noise_cov

def run_lcmv_and_save(subj, task, run, preproc_data_path, source_dir, fwd):
    """
    Apply LCMV beamforming to preprocessed MEG epochs and save epoch-wise
    hippocampal source time series as a .mat file for MATLAB processing.
 
    Two .mat files are saved:
      - *_desc-lcmv-stc_*.mat      : hippocampal source time series per epoch
      - *_desc-lcmv-metadata_*.mat : timing, sampling rate, vertex counts
    """
    # Load preprocessed epochs
    epochs_path = os.path.join(
        preproc_data_path, f"sub-{subj}", 'epochs',
        get_bids_filename(subj, task, run, desc=EPOCHS_DESC)
    )
    print(f"[LCMV] Loading epochs from {epochs_path}")
    epochs_meg = mne.read_epochs(epochs_path, preload=True)
 
    # Check if output already exists
    stc_filename = get_bids_filename(
        subj, task, run, desc='lcmv-stc', suffix='epo', ext='.mat'
    )
    stc_out_path = os.path.join(source_dir, stc_filename)
 
    if os.path.exists(stc_out_path):
        print(f"[LCMV] Output already exists at {stc_out_path}, skipping.")
        return
 
    # Noise covariance — use empty room recording if available, ad hoc otherwise
    print("[LCMV] Computing noise covariance ...")
    noise_cov = compute_noise_cov(epochs_meg, EMPTY_ROOM_PATH)
 
    # Data covariance
    data_cov = compute_covariance(
        epochs_meg, method='shrunk', rank='info',
        n_jobs=LCMV_N_JOBS, verbose=True
    )
 
    # LCMV filter
    print("[LCMV] Computing LCMV filter ...")
    filters = make_lcmv(
        epochs_meg.info, fwd, data_cov,
        reg=LCMV_REG,
        noise_cov=noise_cov,
        reduce_rank=True,
        weight_norm='unit-noise-gain',
        pick_ori='max-power'
    )
 
    # Apply filter to each epoch
    stc_all = apply_lcmv_epochs(epochs_meg, filters)
 
    # Extract hippocampal sources only (skip cortical surface vertices)
    print("[LCMV] Extracting hippocampal source time series ...")
    data_dict = {}
    for idx, stc in enumerate(stc_all):
        n_surf_vertices = len(stc.vertices[0]) + len(stc.vertices[1])
        # Hippocampal sources follow the cortical ones
        hipp_data = stc.data[n_surf_vertices:, :]
        data_dict[f'stc_epoch_{idx}'] = hipp_data
 
    # Save source time series
    print(f"[LCMV] Saving source time series to {stc_out_path}")
    savemat(stc_out_path, {'epochs': data_dict})
 
    # Save metadata (timing, sampling rate, vertex counts per region)
    metadata = {
        'time':          epochs_meg.times,
        'sampling_rate': epochs_meg._raw_sfreq
    }
    for i, src in enumerate(fwd['src']):
        if 'seg_name' in src:
            name = src['seg_name'][:-12] + '_Hippocampus'
        else:
            name = {0: 'Left_Surf', 1: 'Right_Surf'}.get(i, f'source_space_{i}')
        metadata[f'vertno_{name}'] = len(src['vertno'])
 
    meta_filename = get_bids_filename(
        subj, task, run, desc='lcmv-metadata', suffix='epo', ext='.mat'
    )
    meta_out_path = os.path.join(source_dir, meta_filename)
    savemat(meta_out_path, {'metadata': metadata})
    print(f"[LCMV] Metadata saved to {meta_out_path}")

    
#%% MAIN
if __name__ == '__main__':
 
    bids_subj  = f"sub-{SUBJ}"
    source_dir = get_source_dir(SOURCE_RESULTS_PATH, SUBJ)
 
    run_str      = f"_run-{RUN}" if RUN else ""
    raw_fif_path = os.path.join(
        RAW_DATA_PATH,
        f"{bids_subj}_task-{TASK}{run_str}_meg.fif"
    )
 
    print(f"  Hippocampal Source Reconstruction - {bids_subj} - task-{TASK}")
 
    # Step 1: BEM
    bem, mri_subj = compute_or_load_bem(SUBJ, MRI_DIR, source_dir)
 
    # Step 2: Coregistration
    trans_path = load_or_run_coregistration(
        SUBJ, mri_subj, MRI_DIR, source_dir, raw_fif_path
    )
 
    # Step 3: Source space
    mixed_src = compute_or_load_source_space(
        SUBJ, mri_subj, MRI_DIR, source_dir, bem
    )
 
    # Step 4: Forward model
    fwd = compute_or_load_forward(
        SUBJ, mri_subj, source_dir, raw_fif_path,
        trans_path, mixed_src, bem
    )
 
    # Step 5: LCMV + save
    run_lcmv_and_save(
        SUBJ, TASK, RUN, PREPROC_DATA_PATH, source_dir, fwd
    )
 
    print(f"\n[DONE] Results saved to {source_dir}\n")
 