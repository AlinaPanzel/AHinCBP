# Auditory Hyperresponsivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy

[![DOI](https://img.shields.io/badge/OSF-Data-blue)](https://osf.io/pngkr/)
[![Clinical Trial](https://img.shields.io/badge/ClinicalTrials.gov-NCT03294148-green)](https://clinicaltrials.gov/ct2/show/NCT03294148)

**Authors:** Alina E. C. Panzel, Christian Büchel, Andrew Leroux, Tor D. Wager, Yoni K. Ashar

**Language:** MATLAB | **Last updated:** February 2026

---

## Overview

This repository contains all analysis code for a study investigating **auditory hyperresponsivity in chronic back pain (CBP)** and its modifiability through **Pain Reprocessing Therapy (PRT)**.

We compared behavioral and neural responses to aversive auditory and mechanical pressure stimulation in 142 adults with CBP and 51 pain-free controls. CBP patients then entered a randomized controlled trial of PRT vs. open-label placebo vs. usual care ([Ashar et al., 2022, JAMA Psychiatry](https://doi.org/10.1001/jamapsychiatry.2021.2669)).

Analyses include:
- **Behavioral:** Unpleasantness ratings analyzed with linear mixed-effects models (LMMs)
- **Univariate fMRI:** Region-of-interest (ROI) analyses across primary sensory, sensory-integrative (insula), and default mode network regions
- **Multivariate fMRI:** Pattern expression analyses using a priori negative affect patterns ([Čeko et al., 2022](https://doi.org/10.1038/s41593-022-01082-w)) and fibromyalgia-derived patterns ([López-Solà et al., 2017](https://doi.org/10.1097/j.pain.0000000000000707))
- **Longitudinal:** Treatment effects of PRT vs. placebo and usual care on auditory unpleasantness and neural responses
- **Exploratory:** Whole-brain voxel-wise group comparisons and brain-behavior correlations

## Repository Structure

```
AHinCBP/
├── AHinCBP.m                    # Main analysis script (run this)
├── README.md
├── .gitignore
│
├── helper/                      # Helper functions called by AHinCBP.m
│   ├── load_data_objects.m      #   Load & organize fMRI data by group/session
│   ├── load_ROI.m               #   Extract ROI averages from fMRI data
│   ├── load_MVPA.m              #   Apply multivariate pattern classifiers
│   ├── get_behavioral_baseline.m       #   Build baseline behavioral table (long format)
│   ├── get_behavioral_longitudinal.m   #   Build longitudinal behavioral table
│   ├── get_neural_baseline.m           #   Build baseline neural table (ROI + MVPA)
│   ├── get_neural_longitudinal.m       #   Build longitudinal neural table
│   ├── clean_baseline_outliers.m       #   MAD-based outlier removal per cell
│   ├── print_effectsizes.m             #   Compute & print Hedges' g (CBP vs. HC)
│   ├── print_roi_results.m             #   Format ROI result tables
│   ├── pstars.m                        #   Convert p-values to significance stars
│   ├── plot_behavioral_baseline.m      #   Baseline behavioral rating plots
│   ├── plot_baseline_roi.m             #   Baseline ROI activation plots
│   ├── plot_mpfc_treatmenteffects_sound.m            #   mPFC longitudinal plots
│   ├── plot_auditory_treatmenteffects_collapsed.m    #   Treatment effects (collapsed)
│   └── plot_auditory_treatmenteffects_byIntensity.m  #   Treatment effects (by intensity)
│
├── data/                        # Preprocessed data tables (.mat)
│   ├── behavioral_baseline.mat
│   ├── behavioral_longitudinal.mat
│   ├── neural_baseline.mat
│   └── neural_longitudinal_sound.mat
│
└── figures/                     # Generated figures (.png)
    ├── Auditory_Treatmenteffects_collapsed.png
    └── Auditory_Treatmenteffects_byintensity.png
```

## Analysis Pipeline

The full analysis is run sequentially from `AHinCBP.m`, which calls the helper functions in this order:

### Step 1: Load Data
```matlab
d  = get_metadata(d);
[d, lo] = load_data_objects(d, datadir, lo);
a  = load_MSS_atlases(a);
lo = load_ROI(a, lo);
lo = load_MVPA(lo);
```
Loads fMRI contrast images, demographic metadata, ROI atlases, extracts ROI averages, and computes MVPA pattern expression values for all groups and sessions.

### Step 2: Behavioral Analysis (Cross-Sectional)
```matlab
behavioral_baseline = get_behavioral_baseline(d);
lme = fitlme(T, 'Rating ~ cAge + nGender + Intensity*Group + (1|ID)', 'FitMethod','REML');
```
Tests group differences (CBP vs. controls) in unpleasantness ratings for auditory and pressure stimuli using LMMs with Satterthwaite degrees of freedom. Reports main effects of group, intensity, and their interaction.

### Step 3: Neural Baseline Analysis
```matlab
neural_baseline = get_neural_baseline(d, lo);
[neural_baseline_clean, outlier_summary] = clean_baseline_outliers(neural_baseline);
```

- **ROI analyses:** LMMs testing group differences in 11 ROIs (A1, IC, MGN, S1, dIns, vIns, pIns, mPFC, PCC, precuneus) for sound and pressure conditions
- **MVPA analyses:** LMMs testing group differences in 5 multivariate patterns (general negative affect, sound-specific, mechanical-specific, FM-PAIN, FM-MSS)
- **Effect sizes:** Hedges' g for CBP vs. controls at each intensity level
- **Outlier handling:** Median absolute deviation (MAD) removal per ROI/group/modality/intensity cell

### Step 4: Whole-Brain Voxel-Wise Analysis
```matlab
out = regress(dat_mask, 'variable_names', varnames, alpha, 'unc', 'k', kmin, 'robust');
```
Robust regression with weighted group contrast (CBP > HC) on smoothed, gray-matter-masked fMRI data for all four stimulus conditions. Generates t-maps, montages, and cluster tables.

### Step 5: Longitudinal Treatment Effects
```matlab
lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)');
```
Tests group (PRT, placebo, usual care) x timepoint (pre, post) interactions for behavioral ratings and each neural measure, with FDR correction (Benjamini-Hochberg).

## Key Findings

### Cross-Sectional (CBP vs. Controls)

| Domain | Finding | Effect Size |
|--------|---------|-------------|
| **Behavioral** | Greater auditory unpleasantness in CBP | g = 0.95-1.03 |
| **Behavioral** | Greater pressure unpleasantness in CBP | g = 0.49-0.66 |
| **ROI** | Hyperactivation of A1 to sound | g = 0.43-0.66 |
| **ROI** | Hyperactivation of insula (dIns, pIns, vIns) to sound | g = 0.32-0.54 |
| **ROI** | Hypoactivation of mPFC and precuneus to sound | g = -0.22 to -0.58 |
| **MVPA** | Increased general and sound-specific negative affect patterns | g = 0.30-0.39 |
| **MVPA** | Increased fibromyalgia pattern expression (FM-PAIN, FM-MSS) | g = 0.48-0.52 |

### Longitudinal (PRT vs. Controls)

- PRT reduced auditory unpleasantness vs. placebo for low-intensity sounds
- PRT increased mPFC responses vs. usual care, suggesting partial normalization of baseline hypoactivation

## Regions of Interest

ROIs were selected a priori to test three levels of the sensory processing hierarchy:

1. **Primary sensory pathways**
   - Auditory: Primary auditory cortex (A1), inferior colliculus (IC), medial geniculate nucleus (MGN)
   - Somatosensory: Primary somatosensory cortex (S1)

2. **Sensory-integrative pathways**
   - Dorsal-anterior insula (dIns), ventral-anterior insula (vIns), posterior insula (pIns)

3. **Midline default mode network**
   - Medial prefrontal cortex (mPFC), posterior cingulate cortex (PCC), precuneus

**Atlases used:** Chang et al. (2013) for insula; Glasser et al. (2016) for cortical areas; Shen et al. (2013) for IC; Krauth et al. (2010) for MGN.

## Getting Started

### Prerequisites

- **MATLAB R2021b** or later
- **Statistics and Machine Learning Toolbox**
- **Optimization Toolbox**

### External Dependencies

| Toolbox | Purpose | Link |
|---------|---------|------|
| CANlab Core Tools | fMRI data handling, ROI extraction, regression | [canlab.github.io](https://canlab.github.io/repositories/) |
| Measures of Effect Size | Hedges' g, Cohen's d | [MATLAB File Exchange #32398](https://de.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox) |
| al_goodplot | Publication-quality violin/box/scatter plots | [MATLAB File Exchange #91790](https://de.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot) |

### Installation & Usage

1. **Download the data** from the [OSF repository](https://osf.io/pngkr/)

2. **Clone this repository:**
   ```bash
   git clone https://github.com/AlinaPanzel/AHinCBP.git
   ```

3. **Install dependencies** (add to MATLAB path):
   ```matlab
   addpath(genpath('/path/to/CanlabCore'))
   addpath(genpath('/path/to/MeasuresOfEffectSize'))
   addpath(genpath('/path/to/al_goodplot'))
   ```

4. **Configure paths** in `AHinCBP.m`:
   ```matlab
   addpath(genpath('/your/path/to/project'))
   datadir = '/your/path/to/fmri/contrast/maps';
   ```

5. **Run the analysis** step by step in `AHinCBP.m`

### Outputs

- **Console:** LMM summaries, t-statistics, p-values, Hedges' g effect sizes, cluster tables
- **Figures:** Saved to `figures/` (behavioral plots, treatment effect plots)
- **Data tables:** Optionally saved to `data/` as `.mat` files (uncomment save lines in script)

## Statistical Approach

| Analysis | Model | Random Effects |
|----------|-------|----------------|
| Behavioral baseline | `Rating ~ Age + Gender + Intensity * Group + (1\|ID)` | Subject intercept |
| Neural baseline (ROI/MVPA) | `value ~ group + intensity + (1\|subID)` | Subject intercept |
| Longitudinal | `value ~ group * timepoint + intensity + (1 + timepoint\|subID)` | Subject intercept + slope |
| Whole-brain | Robust regression with weighted group contrast | N/A (voxel-wise) |

- Degrees of freedom: Satterthwaite approximation
- Outlier removal: `rmoutliers` (behavioral) and MAD-based (neural)
- Multiple comparisons: FDR (Benjamini-Hochberg) for longitudinal analyses
- Effect sizes: Hedges' g for all pairwise group comparisons

## Related Publication

Panzel, A. E. C., Büchel, C., Leroux, A., Wager, T. D., & Ashar, Y. K. Auditory Hyperresponsivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy. *Annals of Neurology*.

**Clinical Trial:** [NCT03294148](https://clinicaltrials.gov/ct2/show/NCT03294148)

## Data Availability

The data supporting the findings of this study are openly available at: [https://osf.io/pngkr/](https://osf.io/pngkr/)

## Contact

For questions or issues, please open a [GitHub issue](https://github.com/AlinaPanzel/AHinCBP/issues) or contact Alina Panzel.

## Acknowledgements

This work builds on the [CANlab toolbox](https://canlab.github.io/) and source code contributions by Yoni Ashar. Supported by NIH grants R01 DA035484 (NIDA), R01 MH076136 (NIMH), and TL1-TR-002386 (NCATS).

## License

Please contact the authors for licensing information.
