# Auditory Hypersensitivity in Chronic Back Pain:
## A Randomized Controlled Trial of Pain Reprocessing Therapy

Author: Alina Panzel

Last edited: 18.09.2025

Language: MATLAB

## Overview

This repository contains the behavioral and fMRI analysis code for a randomized controlled trial investigating auditory hypersensitivity in chronic back pain (CBP) and the impact of Pain Reprocessing Therapy (PRT).

The project examines behavioral unpleasantness ratings and neural responses to low- and high-intensity aversive auditory and mechanical (thumb pressure) stimulation. Analyses include univariate ROI, multivariate pattern expression (MVPA), and longitudinal treatment effects.

## Analysis Plan

### I. Behavioral Analysis

Unpleasantness ratings for auditory and pressure stimuli.

Linear mixed-effects models test group differences (CBP vs. controls) and group × intensity interactions.

Longitudinal models assess treatment effects (PRT vs. placebo/usual care).

### II. ROI Analysis

Primary sensory regions: Auditory cortex (A1, IC, MGN), Somatosensory cortex (S1).

Sensory-integrative regions: Insula (ventral/dorsal anterior, posterior).

Midline default mode regions: Medial prefrontal cortex (mPFC), Precuneus, Posterior cingulate cortex (PCC).

Group effects tested via LMM

### III. Multivariate Pattern Analysis (MVPA)

Negative affect processing patterns: Generalised (Common), Modality specific (Mechanical, Auditory).

Fibromyalgia-derived patterns: FM_PAIN, FM_MSS.

Group effects tested via LMM ('value ~ group + intensity + (1|subID)')

### IV. Brain–Behavior Correlations

Links neural differences to behavioral sensitivity measures.

### V. Longitudinal Analysis

Tests treatment (PRT vs. placebo/usual care) × timepoint effects.

Treatment effects tested via LMM ('value ~ group*timepoint + intensity + (1 + timepoint | subID)')

## Key Findings 

CBP vs. controls:

Stronger unpleasantness to auditory (g ≈ 0.95–1.03) and mechanical pressure (g ≈ 0.49–0.66).

Auditory fMRI: Hyperactivation in A1 and insula; hypoactivation in mPFC and precuneus.

MVPA: Increased expression of general, auditory, and FM multisensory patterns.

PRT longitudinal effects:

Reduced unpleasantness to low-intensity auditory stimulation.

Increased mPFC activity compared to usual care.

## Usage

Get data at https://neurovault.org/collections/19628/

Clone the repo and open in MATLAB.

Update the addpath and datadir variables in the main script.

Run the analysis step-by-step in AHinCBP.m.

Outputs:

Console: LME model summaries, effect sizes.

Figures: Saved in figures/.

(Optional) Data tables: Saved in data/ if uncommented.

## Dependencies:

MATLAB R2021b+

Statistics & Machine Learning Toolbox

CANlab tool box:         https://canlab.github.io/repositories/

Measure of effect size:  https://de.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox

Plot function:           https://de.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot

## Contact

For questions or issues, please contact Alina Panzel.

## Acknowledgements

This work builds on the CANlab toolbox and source code contributions by Yoni Ashar.
