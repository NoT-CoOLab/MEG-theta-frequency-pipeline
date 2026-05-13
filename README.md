# Hippocampal Theta Frequency Estimation Pipeline

This repository contains a pipeline for estimating hippocampal theta oscillatory frequency from MEG source-reconstructed data.

The pipeline combines:

1. Hippocampal source reconstruction using LCMV beamforming
2. Generalized eigendecomposition (GED) to enhance theta-band activity
3. Cyclic Homogeneous Oscillations (CHO) analysis to detect genuine oscillations and estimate their centre frequency

---

# Repository structure

```text
code/
    run_hippocampal_source_reconstruction.py
    run_theta_frequency_pipeline.m

    functions/
        estimate_hippocampal_theta_frequency.m
        load_hippocampal_sources.m
        define_analysis_windows.m
        run_ged.m
        apply_ged.m
        build_theta_filter.m
        run_cho.m
        retain_frequency_band.m
        ...

external/
    CHO/
        CHO_v22.m
        fft_bandpass_filtering.m
        fft_spectrum.m
        oof_fitting_v2.m
```

---

# Pipeline overview

```text
MEG data
    ↓
LCMV hippocampal source reconstruction
    ↓    
Theta-band GED
    ↓
CHO oscillation detection
    ↓
Theta centre-frequency estimation
```

---

# Dependencies

## Python

- Python 3.10.8
- MNE-Python 1.9.0
- NumPy
- SciPy

## MATLAB

- MATLAB R2023b was used but earlier versions may work
- Signal Processing Toolbox
- Economic Toolbox

---

# External dependencies

## CHO

This pipeline requires the CHO implementation described in:

Cho et al. (2024). *Novel Cyclic Homogeneous Oscillation Detection Method for High Accuracy and Specific Characterization of Neural Dynamics.*

[https://doi.org/10.1101/2023.10.04.560843](https://doi.org/10.1101/2023.10.04.560843)

The CHO implementation is not distributed with this repository.
Please obtain the required CHO MATLAB files from the original authors/source and place them in:

```text
external/CHO/
```

Please cite the original CHO publication when using this pipeline.

---

# How to use

### 1. Clone this repository (or download the files) 

```bash
git clone https://github.com/NoT-CoOLab/MEG-theta-frequency-pipeline.git
```

### 2. Install dependencies
### 3. Obtain the CHO dependency

Please obtain the required CHO MATLAB files from the original authors/source and place them in:

```text
external/CHO/
```
### 4. Run hippocampal source reconstruction

Run:

```bash
python code/run_hippocampal_source_reconstruction.py
```

This produces hippocampal source-reconstructed `.mat` files for each subject.
Expected files per subject:

```text
sub-XX_desc-lcmv-stc.mat
sub-XX_desc-lcmv-metadata.mat
```

### 5. Run theta frequency estimation

Open MATLAB and run:

```matlab
run('code/run_theta_frequency_pipeline.m')
```

Before running:
- edit paths in `run_theta_frequency_pipeline.m`
- ensure all dependencies are on the MATLAB path
- define the parameters 


### Analysis windows

GED and CHO are applied within analysis windows.

Two modes are supported.

### Use loaded epochs directly

```matlab
cfg.window_mode = 'as_loaded';
```

### Fixed-length windows

```matlab
cfg.window_mode = 'fixed';
cfg.window_size = 1000;
```

---

## Output

The pipeline outputs one `.mat` file per subject containing:

```matlab
peak_frequencies
```

Each cell contains detected theta centre frequencies (Hz) for one analysis window.

---

# Citation

If you use this pipeline, please cite:

```text
Marcantoni, E., Daube, C., Wang, D., Cao, C., Sun, B., Zhan, S., Ince, R.A.A., Parkkonen, L., Palva, S., Bush, D., Hanslmayr S. Non-invasive tracking of hippocampal theta oscillations
DOI: https://doi.org/10.64898/2026.01.13.699218
```

Please also cite:
- the original CHO publication
- relevant GED methodological references

---

# Contact
For remaining questions, issues, or suggestions, please use the GitHub [issue tracker](https://github.com/NoT-CoOLab/MEG-theta-frequency-pipeline/issues).
