%% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


% Author: Alina Panzel
% Last Date of Changes: 11.09.2025

%% Load Data & Atlases

% clean 
clc;
clear;

% addpath
addpath(genpath('/Users/alinapanzel/Desktop/Dartmouth_project'))
datadir = '/Users/alinapanzel/Desktop/Dartmouth_project/acute_contrast_maps/work/ics/data/projects/wagerlab/labdata/projects/OLP4CBP/first_level';

% set structs
d = struct; % for directories/ data
lo = struct; % for fmriobjects
a = struct; % for atlases

% Get metadata table
d = get_metadata(d);

% Load in behvaioral & fmri data
[d, lo] = load_data_objects(d, datadir, lo);

% Load atlases
a = load_MSS_atlases(a);


%% Behavioral Analysis

% -------- Baseline Analysis ----------

% Build baseline table in long format (already has centered age and clean vars)
behavioral_baseline = get_baseline_table(d, lo);

% --- AUDIO (Sound) Ratings ---
% Subset only sound trials (Low/High)
AA = behavioral_baseline(behavioral_baseline.Modality=="Sound", :);

% Linear mixed model:
%   Rating ~ Age (centered) + Gender + Intensity (Low/High) * Group (CBP vs Healthy)
%   Random intercept for subject ID
lme_AA = fitlme(AA, ...
    'Rating ~ cAge + nGender + Intensity*Group + (1|ID)', ...
    'FitMethod','REML');

disp('--- Sound model ---');
disp(lme_AA);

% Contrast 1: Group main effect (CBP vs Healthy)
[pval,F,DF1,DF2] = coefTest(lme_AA, [0 1 0 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Sound: Group effect p=%.4f, F=%.2f (df1=%d, df2=%.2f)\n', pval,F,DF1,DF2);

% Contrast 2: Intensity*Group interaction
[pval,F,DF1,DF2] = coefTest(lme_AA, [0 0 0 1 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Sound: Intensity*Group interaction p=%.4f, F=%.2f (df1=%d, df2=%.2f)\n', pval,F,DF1,DF2);


% --- PRESSURE (Thumb) Ratings ---
% Subset only pressure (thumb) trials (Low/High)
PP = behavioral_baseline(behavioral_baseline.Modality=="Pressure", :);

% Linear mixed model
lme_PP = fitlme(PP, ...
    'Rating ~ cAge + nGender + Intensity*Group + (1|ID)', ...
    'FitMethod','REML');

disp('--- Pressure model ---');
disp(lme_PP);

% Contrast 1: Group main effect (CBP vs Healthy)
[pval,F,DF1,DF2] = coefTest(lme_PP, [0 1 0 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Pressure: Group effect p=%.4f, F=%.2f (df1=%d, df2=%.2f)\n', pval,F,DF1,DF2);

% Contrast 2: Intensity main effect (High vs Low)
[pval,F,DF1,DF2] = coefTest(lme_PP, [0 0 1 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Pressure: Intensity effect p=%.4f, F=%.2f (df1=%d, df2=%.2f)\n', pval,F,DF1,DF2);



%% Neural Analysis

% Load in ROI and MVPA values
lo = load_ROI(a, lo);
lo = load_MVPA(lo);
neural_baseline = get_neural_baseline_table(d, lo);
neural_longitudinal = get_longitudinal_table(d,lo);


% ------- Baseline Analysis -------



% ------- Longitudinal Analysis -------




%% Plots & Visualizations

% -------- Baseline Analysis ----------



% ------- Longitudinal Analysis -------





