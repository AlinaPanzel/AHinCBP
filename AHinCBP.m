% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


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

% Baseline analysis 


%% Neural Analysis

% Load in ROI and MVPA values
lo = load_ROI(a, lo);
lo = load_MVPA(lo);


% Test 
test = test

%% Plots & Visualizations

% -------- Baseline Analysis ----------



% ------- Longitudinal Analysis -------

% Get data tables in long format 
con_table_long = get_longformat_table(d,lo, pat_tabl, prexps);





