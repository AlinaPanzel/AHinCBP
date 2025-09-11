% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


% Author: Alina Panzel
% Last Date of Changes: 10.09.2025

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

% Sorts, loads and gets longitudinal fMRI and behavioral metadata 
[d, lo] = get_all_longitudinal(d, lo, datadir);

% Load Atlases
a = load_MSS_atlases(a);


%% Behavioral Analysis

% Baseline analysis 


%% Neural Analysis


% Load ROIs of comparison interest
areas = {'A1','MGN','IC','S1_L','S1_R', 'm_ventral_insula', 'm_dorsal_insula', 'm_posterior_insula', 'mPFC', 'precuneus', 'PCC'};
sessions = {'S1','S2'};
Treatment = {'G1' 'G2' 'G3'};
field =  {'all_s_l', 'all_s_h', 'completers_s_l', 'completers_s_h'}; %'s_l' 's_h', 
for n = 1:numel(Treatment)
    for s = 1:numel(sessions)
        for u = 1:numel(areas)
            for i = 1:numel(field)  %numel(fields) % number of conditions
                r = extract_roi_averages(lo.(Treatment{n}).(sessions{s}).(field{i}), a.(areas{u})); 
                lo.roi.(Treatment{n}).(sessions{s}).(areas{u}).(field{i}) = cat(2, r.dat); 
            end
        end
    end
end

% Load MVPA patterns of comparison interest
for g = 1:numel(Treatment)
    for s = 1:numel(sessions)
        for i = 1:numel(field)
            pat_tabl.(Treatment{g}).(sessions{s}).(field{i}) = apply_FM_patterns(lo.(Treatment{g}).(sessions{s}).(field{i}), 'cosine_similarity'); % Fibromyalgia patterns
            prexps.(Treatment{g}).(sessions{s}).(field{i}) = apply_multiaversive_mpa2_patterns(lo.(Treatment{g}).(sessions{s}).(field{i})); % Negative affect patterns
        end
    end 
end



% Baseline Analysis 


% Longitudinal Analysis 

% Get data tables in long format 
con_table_long = get_longformat_table(d,lo, pat_tabl, prexps);









% Get sorted brain measures -> only wors
%[con_table_lo, con_table_hi] = get_data_tables(d, pat_tabl, prexps)
