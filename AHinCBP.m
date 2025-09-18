%% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


% Author: Alina Panzel
% Last Date of Changes: 13.09.2025

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

% Build baseline table in long format (has centered age and clean vars)
behavioral_baseline = get_behavioral_baseline(d); 
behavioral_longitudinal = get_behavioral_longitudinal(d);
%save(fullfile('data','behavioral_baseline.mat'), 'behavioral_baseline')
%save(fullfile('data','behavioral_longitudinal.mat'), 'behavioral_longitudinal')


% --- AUDIO (Sound) Ratings ---

% Subset only sound trials (Low/High)
AA = behavioral_baseline(behavioral_baseline.Modality=="Sound", :);

% Remove outliers on Rating
[~, idxOut] = rmoutliers(AA.Rating);
AA(idxOut,:) = [];

%histogram(AA.Rating) --> none removed

% Linear mixed model:
%   Rating ~ Age (centered) + Gender + Intensity (Low/High) * Group (CBP vs Healthy)
%   Random intercept for subject ID
lme_AA = fitlme(AA, 'Rating ~ cAge + nGender + Intensity*Group + (1|ID)','FitMethod','REML');

disp('--- Sound model ---'); disp(lme_AA);

% Group main effect
[p,F,DF1,DF2] = coefTest(lme_AA, [0 0 1 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Sound: Group effect        F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);

% Intensity main effect
[p,F,DF1,DF2] = coefTest(lme_AA, [0 1 0 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Sound: Intensity effect    F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);

% Interaction
[p,F,DF1,DF2] = coefTest(lme_AA, [0 0 0 0 0 1], 0, 'DFMethod','Satterthwaite');
fprintf('Sound: Intensity×Group     F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);



% --- PRESSURE (Thumb) Ratings ---
% Subset only pressure (thumb) trials (Low/High)
PP = behavioral_baseline(behavioral_baseline.Modality=="Pressure", :);

% remove outliers on Rating
[~, idxOut] = rmoutliers(PP.Rating);
PP(idxOut,:) = [];

lme_PP = fitlme(PP, ...
    'Rating ~ cAge + nGender + Intensity*Group + (1|ID)', ...
    'FitMethod','REML');

disp('--- Pressure model ---'); disp(lme_PP);

% Group main effect
[p,F,DF1,DF2] = coefTest(lme_PP, [0 0 1 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Pressure: Group effect     F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);

% Intensity main effect
[p,F,DF1,DF2] = coefTest(lme_PP, [0 1 0 0 0 0], 0, 'DFMethod','Satterthwaite');
fprintf('Pressure: Intensity effect F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);

% Interaction
[p,F,DF1,DF2] = coefTest(lme_PP, [0 0 0 0 0 1], 0, 'DFMethod','Satterthwaite');
fprintf('Pressure: Intensity×Group  F(%d, %.2f)=%.2f, p=%.4g\n', DF1, DF2, F, p);



% -------- Longitudinal Analysis ----------

% PLot treatment effects
f1 = plot_auditory_treatmenteffects_collapsed(d);
saveas(f1, fullfile('figures','auditory_treatmenteffects_collapsed.png'));
exportgraphics(f1, fullfile('figures','auditory_treatmenteffects_collapsed.png'), ...
               'Resolution',300);

f2 = plot_auditory_treatmenteffects_byIntensity(d);




%% Neural Analysis

% Load in ROI and MVPA values
lo = load_ROI(a, lo);
lo = load_MVPA(lo);

% Get longformat data tables
neural_baseline = get_neural_baseline(d, lo);
neural_longitudinal_sound = get_neural_longitudinal(d,lo);

% Save data
%save(fullfile('data','neural_baseline.mat'), 'neural_baseline')
%save(fullfile('data','neural_longitudinal_sound.mat'), 'neural_longitudinal')

% DO outlier removal beforehand
[neural_baseline_clean, outlier_summary] = clean_baseline_outliers(neural_baseline);


%% %% Neural Baseline: LMM SOUND

fprintf('\n==== LMM RESULTS SOUND ====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

res_sound = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','df1','df2','p_group'});

dropmvpa = {'general','sound','FM_PAIN','FM_MSS'};

neural_baseline_roi_sound = neural_baseline_clean( ...
    neural_baseline_clean.modality == "Sound" & ...
    ~ismember(neural_baseline_clean.measure, dropmvpa), :);

allROIs = cellstr(unique(neural_baseline_roi_sound.measure));

for i = 1:numel(allROIs)
    R  = allROIs{i};

    Tk = neural_baseline_roi_sound(neural_baseline_roi_sound.measure==R, :);

    % group×intensity means (just for printing)
    G = groupsummary(Tk, {'GroupBin','intensity'}, 'mean', 'value');
    % make sure we have all cells; fill missing with NaN
    gcat = {'HC','CBP'}; icat = {'Low','High'};
    getMean = @(g,i) mean(G.mean_value(G.GroupBin==g & G.intensity==i),'omitnan');

    hcL = getMean('HC','Low');  hcH = getMean('HC','High');
    cbL = getMean('CBP','Low'); cbH = getMean('CBP','High');

    % n per group (unique subjects after cleaning)
    nHC  = numel(unique(Tk.subID(Tk.GroupBin=="HC"))); % for double check/oprinting
    nCBP = numel(unique(Tk.subID(Tk.GroupBin=="CBP")));

    % --- LMM ---
    lme = fitlme(Tk, 'value ~ GroupBin + intensity + (1|subID)', 'FitMethod','REML');
    [~,~,stats] = fixedEffects(lme, 'DFMethod','Satterthwaite'); % the DF are inflated, this corrects that
    
    % the effect of interest is Group. This is the 3rd row.
    p = stats.pValue(3);
    DF1 = stats.DF(3);
    DF2 = NaN; % why should there be two DFs? 
    T = stats.tStat(3); % T not F; not renaming now

    % print one-line summary ROI
    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        R, nHC, nCBP, hcL, hcH, cbL, cbH, T, DF1, p, pstars(p)); % pvals = Group effect


    res_sound = [res_sound; {string(R), height(Tk), T, DF1, DF2, p}]; 
end


%% Neural Baseline: LMM PRESSURE

fprintf('\n==== LMM RESULTS PRESSURE====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

dropmvpa = {'general','sound','FM_PAIN','FM_MSS'};

res_pressure = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','df1','df2','p_group'});

neural_baseline_roi_pressure = neural_baseline_clean( ...
    neural_baseline_clean.modality == "Pressure" & ...
    ~ismember(neural_baseline_clean.measure, dropmvpa), :);

for i = 1:numel(allROIs)
    R  = allROIs{i};

    Tk = neural_baseline_roi_pressure(neural_baseline_roi_pressure.measure==R, :);
    
    % group×intensity means (for printing)
    G = groupsummary(Tk, {'GroupBin','intensity'}, 'mean', 'value');
    % ensure we have all cells; fill missing with NaN
    gcat = {'HC','CBP'}; icat = {'Low','High'};
    getMean = @(g,i) mean(G.mean_value(G.GroupBin==g & G.intensity==i),'omitnan');

    hcL = getMean('HC','Low');  hcH = getMean('HC','High');
    cbL = getMean('CBP','Low'); cbH = getMean('CBP','High');

    % n per group (unique subjects after cleaning)
    nHC  = numel(unique(Tk.subID(Tk.GroupBin=="HC")));
    nCBP = numel(unique(Tk.subID(Tk.GroupBin=="CBP")));

    % --- LMM ---
    lme = fitlme(Tk, 'value ~ GroupBin + intensity + (1|subID)', 'FitMethod','REML');
    [~,~,stats] = fixedEffects(lme, 'DFMethod','Satterthwaite'); % the DF are inflated, this corrects that
    
    % the effect of interest is Group. This is the 3rd row.
    p = stats.pValue(3);
    DF1 = stats.DF(3);
    DF2 = NaN; % why should there be two DFs? 
    T = stats.tStat(3); % T not F; not renaming now

    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        R, nHC, nCBP, hcL, hcH, cbL, cbH, T, DF1, p, pstars(p));

    res_pressure = [res_pressure; {string(R), height(Tk), T, DF1, DF2, p}]; 
end


%% Neural Longitudinal: LMM SOUND 

% Fit & test treatment effects with LMM (random slopes & intercepts)

% Set up & clearn
T = neural_longitudinal_sound;                             % columns: subID,timepoint,group,intensity,measure,value
T.value     = double(T.value);
T.subID     = categorical(T.subID);
T.group     = categorical(T.group);
T.timepoint = categorical(T.timepoint);
T.intensity = categorical(T.intensity);
T.measure   = categorical(T.measure);

measures = categories(T.measure);

results = table(strings(0,1), zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'measure','NumObs','F','df1','df2','p'});
models  = struct('measure',{},'lme',{});

for k = 1:numel(measures)
    m  = measures{k};
    Tk = T(T.measure == m, :);
    
    % Make sure predictors are categorical/numeric as needed
    Tk.group     = categorical(Tk.group);
    Tk.timepoint = double(Tk.timepoint);   % 1=baseline, 2=post
    Tk.intensity = double(Tk.intensity);   % 1=low, 2=high
    
    % Random intercept + random slope for timepoint
    lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)');
    
    % Joint test of group × timepoint
    cn = string(lme.CoefficientNames);
    J  = find(contains(cn,"group") & contains(cn,"timepoint"));
    
    if ~isempty(J)
        L = zeros(numel(J), numel(cn));
        for i=1:numel(J), L(i,J(i)) = 1; end
        [p,F,df1,df2] = coefTest(lme, L);
    else
        F=NaN; df1=NaN; df2=NaN; p=NaN;
    end

    results = [results; {string(m), height(Tk), F, df1, df2, p}];
    models(end+1).measure = string(m); %#ok<SAGROW>
    models(end).lme = lme;
end

% FDR (Benjamini–Hochberg)
p = results.p; keep = ~isnan(p); q = nan(size(p));
if any(keep), [~,~,~,q(keep)] = fdr_bh(p(keep)); end
results.q_BH = q;

% nice printing table
Rprint = sortrows(results, 'p');
Rprint.F   = round(Rprint.F,3);
Rprint.df1 = round(Rprint.df1,2);
Rprint.df2 = round(Rprint.df2,2);
Rprint.p   = round(Rprint.p,4);
Rprint.q_BH= round(Rprint.q_BH,4);

disp('=== Group × Time interaction (per measure) ===');
disp(Rprint);

