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
behavioral_baseline = get_baseline_table(d, lo);
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

neural_baseline_roi = get_neural_baseline_roi(d, lo);
neural_baseline_mvpa = get_neural_baseline_mvpa(d,lo);
neural_longitudinal = get_longitudinal_table(d,lo);

% Save data
%save(fullfile('data','neural_baseline.mat'), 'neural_baseline')
%save(fullfile('data','neural_longitudinal.mat'), 'neural_longitudinal')


%% Neural Baseline: Independent T-test

% SOUND
sound_baseline_roi = neural_baseline_roi(neural_baseline_roi.modality == "Sound",:); 

T = sound_baseline_roi;
T.GroupBin = categorical(ismember(T.group,[1 2 3]),[0 1],{'HC','CBP'});

ROIs = categories(categorical(T.measure));
ints = categories(categorical(T.intensity));

for ii = 1:numel(ints)
    thisInt = ints{ii};
    fprintf('\n=== Intensity: %s ===\n', string(thisInt));

    p = nan(numel(ROIs),1);
    statsOut = cell(numel(ROIs),1);

    for r = 1:numel(ROIs)
        Tk = T(T.measure==ROIs{r} & T.intensity==thisInt,:);
        x = Tk.value(Tk.GroupBin=="HC");
        y = Tk.value(Tk.GroupBin=="CBP");
        x = rmoutliers(x);
        y = rmoutliers(y);
        [~,p(r),~,st] = ttest2(x,y);
        statsOut{r} = sprintf('%s: nHC=%d, nCBP=%d, meanHC=%.2f, meanCBP=%.2f, t=%.2f, p=%.4f',...
            ROIs{r},numel(x),numel(y),mean(x,'omitnan'),mean(y,'omitnan'),st.tstat,p(r));
    end

    % FDR correction
   [~,~,~,q] = fdr_bh(p); % if you have fdr_bh
   
    % Print results
    for r=1:numel(ROIs)
        fprintf('%s, q=%.4f\n', statsOut{r}, q(r));
    end
end



% PRESSURE

pressure_baseline_roi = neural_baseline_roi(neural_baseline_roi.modality == "Pressure",:); 

T = pressure_baseline_roi;

ROIs = categories(categorical(T.measure));
ints = categories(categorical(T.intensity));

for ii = 1:numel(ints)
    thisInt = ints{ii};
    fprintf('\n=== Intensity: %s ===\n', string(thisInt));

    p = nan(numel(ROIs),1);
    statsOut = cell(numel(ROIs),1);

    for r = 1:numel(ROIs)
        Tk = T(T.measure==ROIs{r} & T.intensity==thisInt,:);
        x = Tk.value(Tk.GroupBin=="HC");
        y = Tk.value(Tk.GroupBin=="CBP");
        x = rmoutliers(x);
        y = rmoutliers(y);
        [~,p(r),~,st] = ttest2(x,y);
        statsOut{r} = sprintf('%s: nHC=%d, nCBP=%d, meanHC=%.2f, meanCBP=%.2f, t=%.2f, p=%.4f',...
            ROIs{r},numel(x),numel(y),mean(x,'omitnan'),mean(y,'omitnan'),st.tstat,p(r));
    end

    % FDR correction
   [~,~,~,q] = fdr_bh(p); % if you have fdr_bh
   
    % Print results
    for r=1:numel(ROIs)
        fprintf('%s, q=%.4f\n', statsOut{r}, q(r));
    end
end


%% Neural Baseline: LMM

% Binary group for modelling/printing: HC vs CBP
neural_baseline_roi.GroupBin = categorical(ismember(neural_baseline_roi.group,[1 2 3]), [0 1], {'HC','CBP'});

% Ensure intensity is categorical Low/High
if ~iscategorical(neural_baseline_roi.intensity)
    neural_baseline_roi.intensity = categorical(neural_baseline_roi.intensity,[1 2],{'Low','High'});
end

allROIs = categories(categorical(neural_baseline_roi.measure));

% SOUND
fprintf('\n==== SOUND analysis (all ROIs) ====\n');
res_sound = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','F_group','df1','df2','p_group'});

for i = 1:numel(allROIs)
    R  = allROIs{i};

    % --- subset this ROI for SOUND only (keep if you’re separating modalities this way) ---
    Tk = neural_baseline_roi(neural_baseline_roi.measure==R & ...
                             neural_baseline_roi.modality=="Sound", :);  % if you have a 'modality' var
    if isempty(Tk), continue; end

    % --- per-group outlier removal on 'value' ---
    Lvls = categories(Tk.GroupBin);
    idxDrop = false(height(Tk),1);
    removedCounts = zeros(numel(Lvls),1);

    for g = 1:numel(Lvls)
        gi = Tk.GroupBin == Lvls{g};
        if any(gi)
            [~, out] = rmoutliers(Tk.value(gi));
            removedCounts(g) = sum(out);
            if removedCounts(g)>0
                % mark rows to drop (only those flagged within this group)
                idxLocal = find(gi);
                idxDrop(idxLocal(out)) = true;
            end
        end
    end
    Tk(idxDrop,:) = [];

    if height(Tk) < 4, continue; end

    % --- print group means (by intensity) and outlier counts ---
    Gstats = groupsummary(Tk, {'GroupBin','intensity'}, 'mean', 'value');
    fprintf('\nROI: %s (SOUND)\n', R);
    disp(Gstats(:, {'GroupBin','intensity','mean_value'}));
    for g = 1:numel(Lvls)
        fprintf('  Outliers removed in %-4s: %d\n', char(Lvls{g}), removedCounts(g));
    end

    % --- LMM (group main effect) ---
    lme = fitlme(Tk, 'value ~ GroupBin*intensity + (1|subID)', 'FitMethod','REML');
    [p,F,DF1,DF2] = coefTest(lme, [0 1 0 0], 0, 'DFMethod','Satterthwaite');
    fprintf('  Group effect: F(%d,%.2f)=%.2f, p=%.4g   (N=%d rows)\n', DF1, DF2, F, p, height(Tk));

    res_sound = [res_sound; {string(R), height(Tk), F, DF1, DF2, p}]; 
end


% PRESSURE
fprintf('\n==== PRESSURE analysis (all ROIs) ====\n');
res_press = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','F_group','df1','df2','p_group'});

for i = 1:numel(allROIs)
    R  = allROIs{i};

    % --- subset this ROI for PRESSURE only ---
    Tk = neural_baseline_roi(neural_baseline_roi.measure==R & ...
                             neural_baseline_roi.modality=="Pressure", :);  
    if isempty(Tk), continue; end

    % --- per-group outlier removal on 'value' ---
    Lvls = categories(Tk.GroupBin);
    idxDrop = false(height(Tk),1);
    removedCounts = zeros(numel(Lvls),1);

    for g = 1:numel(Lvls)
        gi = Tk.GroupBin == Lvls{g};
        if any(gi)
            [~, out] = rmoutliers(Tk.value(gi));
            removedCounts(g) = sum(out);
            if removedCounts(g)>0
                idxLocal = find(gi);
                idxDrop(idxLocal(out)) = true;
            end
        end
    end
    Tk(idxDrop,:) = [];

    if height(Tk) < 4, continue; end

    % --- print group means (by intensity) and outlier counts ---
    Gstats = groupsummary(Tk, {'GroupBin','intensity'}, 'mean', 'value');
    fprintf('\nROI: %s (PRESSURE)\n', R);
    disp(Gstats(:, {'GroupBin','intensity','mean_value'}));
    for g = 1:numel(Lvls)
        fprintf('  Outliers removed in %-4s: %d\n', char(Lvls{g}), removedCounts(g));
    end

    % --- LMM (group main effect) ---
    lme = fitlme(Tk, 'value ~ GroupBin*intensity + (1|subID)', 'FitMethod','REML');
    [p,F,DF1,DF2] = coefTest(lme, [0 1 0 0], 0, 'DFMethod','Satterthwaite');
    fprintf('  Group effect: F(%d,%.2f)=%.2f, p=%.4g   (N=%d rows)\n', DF1, DF2, F, p, height(Tk));

    res_press = [res_press; {string(R), height(Tk), F, DF1, DF2, p}]; %#ok<AGROW>
end




%% Neural Longitudinal 

% Fit & test treatment effects with LMM (random slopes & intercepts)

% Set up & clearn
T = neural_longitudinal;                             % columns: subID,timepoint,group,intensity,measure,value
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


% FOR YONI TO CHECK

Tk = T(T.measure == 'mPFC', :);
lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)');
lme

Tk = T(T.measure == 'precuneus', :);
lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)');
lme

% Not fully functional yet
%visualize_lme_interaction(T, 'mPFC');
%visualize_lme_interaction(T, 'precuneus');

%plot_treatmenteffects();



%% Plots & Visualizations

% -------- Baseline Analysis ----------



% ------- Longitudinal Analysis -------





