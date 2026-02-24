%% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


% Author: Alina Panzel
% Last Date of Changes: 23.02.2026

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
AA = behavioral_baseline(behavioral_baseline.Modality== "Sound", :);

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


% -------- Plot Unpleasantness -----------

% Define y-limits for the plots
yLimit = [-10 110];

% Plot data for different regions (figure created inside plot_figures)
plot_figures(d, [1 2 3 4], 'Unpleasantness Ratings','behavioral', yLimit, 1, false);


% -------- VAS Characteristics ---------

% Sort & retrieve means and variance of spontaneous pain and unpleasantness
% ratings
baseline_ratings = get_rating_characteristics(d);



% -------- Behavioural Correlations ----------

% --- Back Pain (BPI-SF) vs. Auditory Unpleasantness (CBP only) ---

% Subset CBP patients at baseline with auditory ratings
M = d.metadata;
M = M(M.time == 1, :);                           % baseline only
isCBP = ismember(M.group, [1 2 3]);
M_cbp = M(isCBP, :);

fprintf('\n==== BACK PAIN vs. AUDITORY UNPLEASANTNESS (CBP only) ====\n');

% Low intensity sound
valid = ~isnan(M_cbp.pain_avg) & ~isnan(M_cbp.acute_mean_sound_lo);
x = M_cbp.pain_avg(valid);
y = M_cbp.acute_mean_sound_lo(valid);
[r, p] = corr(x, y);
fprintf('Sound Low:  r(%d) = %.2f, p = %.3g%s\n', numel(x)-2, r, p, pstars(p));

% High intensity sound
valid = ~isnan(M_cbp.pain_avg) & ~isnan(M_cbp.acute_mean_sound_hi);
x = M_cbp.pain_avg(valid);
y = M_cbp.acute_mean_sound_hi(valid);
[r, p] = corr(x, y);
fprintf('Sound High: r(%d) = %.2f, p = %.3g%s\n', numel(x)-2, r, p, pstars(p));


% --- Spontaneous Pain vs. Task-Evoked Unpleasantness (CBP only) ---

fprintf('\n==== SPONTANEOUS PAIN vs. TASK-EVOKED UNPLEASANTNESS (CBP only) ====\n');

% Compute spontaneous pain variance from raw ratings where available
if ismember('spon_pain_ratings', M_cbp.Properties.VariableNames)
    spon_var = nan(height(M_cbp),1);
    for si = 1:height(M_cbp)
        ratings_i = M_cbp.spon_pain_ratings(si,:);
        if any(~isnan(ratings_i))
            spon_var(si) = nanvar(ratings_i);
        end
    end
    M_cbp.spon_pain_var = spon_var;
end

% Four conditions: sound_lo, sound_hi, pressure_lo, pressure_hi
condLabels = {'Sound Low','Sound High','Pressure Low','Pressure High'};
condVars   = {'acute_mean_sound_lo','acute_mean_sound_hi', ...
              'acute_mean_thumb_lo','acute_mean_thumb_hi'};

fprintf('\n--- Mean spontaneous pain ---\n');
fprintf('%-16s %8s %10s %6s\n', 'Condition','r','p','');
fprintf('%s\n', repmat('-',1,44));

for ci = 1:numel(condVars)
    valid = ~isnan(M_cbp.spon_pain_rating_mean) & ~isnan(M_cbp.(condVars{ci}));
    x = M_cbp.spon_pain_rating_mean(valid);
    y = M_cbp.(condVars{ci})(valid);
    [r, p] = corr(x, y);
    fprintf('%-16s %8.2f %10.3g %6s\n', condLabels{ci}, r, p, pstars(p));
end

if ismember('spon_pain_var', M_cbp.Properties.VariableNames)
    fprintf('\n--- Variance of spontaneous pain ---\n');
    fprintf('%-16s %8s %10s %6s\n', 'Condition','r','p','');
    fprintf('%s\n', repmat('-',1,44));

    for ci = 1:numel(condVars)
        valid = ~isnan(M_cbp.spon_pain_var) & ~isnan(M_cbp.(condVars{ci}));
        x = M_cbp.spon_pain_var(valid);
        y = M_cbp.(condVars{ci})(valid);
        [r, p] = corr(x, y);
        fprintf('%-16s %8.2f %10.3g %6s\n', condLabels{ci}, r, p, pstars(p));
    end
end



% -------- Behavioural Longitudinal Analysis ----------

% PLot treatment effects
f1 = plot_auditory_treatmenteffects_collapsed(d);
saveas(f1, fullfile('figures','auditory_treatmenteffects_collapsed.png'));
exportgraphics(f1, fullfile('figures','auditory_treatmenteffects_collapsed.png'), ...
    'Resolution',300);

plot_auditory_treatmenteffects_byIntensity(d);




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

% Do outlier removal beforehand
[neural_baseline_clean, outlier_summary] = clean_baseline_outliers(neural_baseline);

%% Neural Baseline: LMM ROI

%==== SOUND ====

fprintf('\n==== LMM RESULTS SOUND ====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

res_sound = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','DF','p_group'});

dropmvpa = {'general','sound','mechanical','FM_PAIN','FM_MSS'};

neural_baseline_roi_sound = neural_baseline_clean( ...
    neural_baseline_clean.modality == "Sound" & ...
    ~ismember(neural_baseline_clean.measure, dropmvpa), :);

allROIs = cellstr(unique(neural_baseline_roi_sound.measure));

for i = 1:numel(allROIs)
    R  = allROIs{i};

    Tk = neural_baseline_roi_sound(neural_baseline_roi_sound.measure==R, :);

    % group×intensity means (just for printing)
    G = groupsummary(Tk, {'group','intensity'}, 'mean', 'value');
    % make sure we have all cells; fill missing with NaN
    gcat = {'HC','CBP'}; icat = {'Low','High'};
    getMean = @(g,i) mean(G.mean_value(G.group==g & G.intensity==i),'omitnan');

    hcL = getMean('HC','Low');  hcH = getMean('HC','High');
    cbL = getMean('CBP','Low'); cbH = getMean('CBP','High');

    % n per group (unique subjects after cleaning)
    nHC  = numel(unique(Tk.subID(Tk.group=="HC"))); % for double check/oprinting
    nCBP = numel(unique(Tk.subID(Tk.group=="CBP")));

    % --- LMM ---
    lme = fitlme(Tk, 'value ~ group + intensity + (1|subID)', 'FitMethod','REML');
    [~,~,stats] = fixedEffects(lme, 'DFMethod','Satterthwaite'); % the DF are inflated, this corrects that

    % Find the coefficient for CBP vs HC regardless of ordering -> did this
    % because I restructured the tables
    ix = find(strcmp(stats.Name,'group'));
    if isempty(ix)
        ix = find(contains(stats.Name,'group'));
    end
    t  = stats.tStat(ix);
    df = stats.DF(ix);
    p  = stats.pValue(ix);

    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        R, nHC, nCBP, hcL, hcH, cbL, cbH, t, df, p, pstars(double(p)));


    res_sound = [res_sound; {string(R), height(Tk), t, df, p}];
end

print_effectsizes(neural_baseline_roi_sound, 'ROI (Sound)');


%===== PRESSURE =====

fprintf('\n==== LMM RESULTS PRESSURE====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

dropmvpa = {'general','sound','mechanical','FM_PAIN','FM_MSS'};

res_pressure = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','DF','p_group'});

neural_baseline_roi_pressure = neural_baseline_clean( ...
    neural_baseline_clean.modality == "Pressure" & ...
    ~ismember(neural_baseline_clean.measure, dropmvpa), :);

for i = 1:numel(allROIs)
    R  = allROIs{i};

    Tk = neural_baseline_roi_pressure(neural_baseline_roi_pressure.measure==R, :);

    % group×intensity means (for printing)
    G = groupsummary(Tk, {'group','intensity'}, 'mean', 'value');
    % ensure we have all cells; fill missing with NaN
    gcat = {'HC','CBP'}; icat = {'Low','High'};
    getMean = @(g,i) mean(G.mean_value(G.group==g & G.intensity==i),'omitnan');

    hcL = getMean('HC','Low');  hcH = getMean('HC','High');
    cbL = getMean('CBP','Low'); cbH = getMean('CBP','High');

    % n per group (unique subjects after cleaning)
    nHC  = numel(unique(Tk.subID(Tk.group=="HC")));
    nCBP = numel(unique(Tk.subID(Tk.group=="CBP")));

    % --- LMM ---
    lme = fitlme(Tk, 'value ~ group + intensity + (1|subID)', 'FitMethod','REML');
    [~,~,stats] = fixedEffects(lme, 'DFMethod','Satterthwaite'); % the DF are inflated, this corrects that

    % Find the coefficient for CBP vs HC regardless of ordering -> did this
    % because I restructured the tables
    ix = find(strcmp(stats.Name,'group'));
    if isempty(ix)
        ix = find(contains(stats.Name,'group'));
    end
    t  = stats.tStat(ix);
    df = stats.DF(ix);
    p  = stats.pValue(ix);
   
    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        R, nHC, nCBP, hcL, hcH, cbL, cbH, t, df, p, pstars(double(p)));
    res_pressure = [res_pressure; {string(R), height(Tk), t, df, p}];
end

print_effectsizes(neural_baseline_roi_pressure, 'ROI (PRESSURE)');

% plot ROI

get_roiplots(lo)

%% Neural Baseline: LMM MVPA

%====== SOUND ======

fprintf('\n==== LMM RESULTS SOUND (MVPAs) ====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'MVPA','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group t [df] (p)');
fprintf('%s\n', repmat('-',1,90));

mvpaKeep = {'general','sound','mechanical','FM_PAIN','FM_MSS'};

% subset + light type coercion
T = neural_baseline_clean( ...
        neural_baseline_clean.modality=="Sound" & ...
        ismember(neural_baseline_clean.measure, mvpaKeep), :);

T.value     = double(T.value);
T.group     = categorical(string(T.group), {'HC','CBP'});     % lock refs: HC baseline
T.intensity = categorical(string(T.intensity), {'Low','High'}); % lock refs: Low baseline
T.measure   = removecats(categorical(string(T.measure)));

allMVPAs = cellstr(categories(T.measure));

% ---- FIX: 5 columns for 5 names ----
res_sound = table( ...
    string.empty, ...   % MVPA
    zeros(0,1), ...     % N
    nan(0,1), ...       % t_group
    nan(0,1), ...       % df
    nan(0,1), ...       % p_group
    'VariableNames', {'MVPA','N','t_group','df','p_group'});

for i = 1:numel(allMVPAs)
    M  = allMVPAs{i};
    Tk = T(T.measure==M, :);
    if isempty(Tk), continue; end

    % means for printing
    G = groupsummary(Tk, {'group','intensity'}, 'mean', 'value');
    getMean = @(g,i) mean(G.mean_value(G.group==g & G.intensity==i),'omitnan');
    hcL = getMean('HC','Low');   hcH = getMean('HC','High');
    cbL = getMean('CBP','Low');  cbH = getMean('CBP','High');

    % n per group (unique subjects)
    nHC  = numel(unique(Tk.subID(Tk.group=="HC")));
    nCBP = numel(unique(Tk.subID(Tk.group=="CBP")));

    % LMM: Group main effect (averaged over intensity), with fixed coding
    lme = fitlme(Tk, 'value ~ group + intensity + (1|subID)', ...
                 'FitMethod','REML','DummyVarCoding','reference');

    % pull the Group coefficient by NAME (stable to ordering)
    [~,~,FE] = fixedEffects(lme,'DFMethod','Satterthwaite');
    ix = find(strcmp(FE.Name,'group_CBP'));     % CBP vs HC (HC is ref)
    if isempty(ix)
        ix = find(contains(FE.Name,'group'));   % fallback if label formatting differs
    end
    t  = FE.tStat(ix);
    df = FE.DF(ix);
    p  = FE.pValue(ix);

    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        M, nHC, nCBP, hcL, hcH, cbL, cbH, t, df, p, pstars(p));

    res_sound = [res_sound; {string(M), height(Tk), t, df, p}]; %#ok<AGROW>
end

T_sound = T;
print_effectsizes(T_sound, 'MVPA (Sound)');


% ===== PRESSURE =====

fprintf('\n==== LMM RESULTS PRESSURE (MVPAs) ====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'MVPA','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group t [df] (p)');
fprintf('%s\n', repmat('-',1,90));

T = neural_baseline_clean( ...
    neural_baseline_clean.modality=="Pressure" & ...
    ismember(neural_baseline_clean.measure, mvpaKeep), :);

T.value     = double(T.value);
T.group     = categorical(string(T.group), {'HC','CBP'});     % lock refs: HC baseline
T.intensity = categorical(string(T.intensity), {'Low','High'}); % lock refs: Low baseline
T.measure   = removecats(categorical(string(T.measure)));

allMVPAs = cellstr(categories(T.measure));

%Preallocate
res_pressure = table(string.empty,zeros(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'MVPA','N','t_group','df','p_group'});

for i = 1:numel(allMVPAs)
    M  = allMVPAs{i};
    Tk = T(T.measure==M, :);
    if isempty(Tk), continue; end

    % means for printing
    G = groupsummary(Tk, {'group','intensity'}, 'mean', 'value');
    getMean = @(g,i) mean(G.mean_value(G.group==g & G.intensity==i),'omitnan');
    hcL = getMean('HC','Low');   hcH = getMean('HC','High');
    cbL = getMean('CBP','Low');  cbH = getMean('CBP','High');

    % n per group (unique subjects)
    nHC  = numel(unique(Tk.subID(Tk.group=="HC")));
    nCBP = numel(unique(Tk.subID(Tk.group=="CBP")));

    % LMM: Group main effect (averaged over intensity), with fixed coding
    lme = fitlme(Tk, 'value ~ group + intensity + (1|subID)', ...
        'FitMethod','REML','DummyVarCoding','reference');

    % pull the Group coefficient by NAME (stable to ordering)
    [~,~,stats] = fixedEffects(lme,'DFMethod','Satterthwaite');
    ix = find(strcmp(stats.Name,'group_CBP'));     % CBP vs HC (HC is ref)
    if isempty(ix)
        ix = find(contains(stats.Name,'group'));   % fallback if label formatting differs
    end
    t  = stats.tStat(ix);
    df = stats.DF(ix);
    p  = stats.pValue(ix);

    fprintf('%-16s | %4d %4d | %8.2f %8.2f | %8.2f %8.2f | %6.2f [%3.2f] (%.3g)%s\n', ...
        M, nHC, nCBP, hcL, hcH, cbL, cbH, t, df, p, pstars(p));

    res_pressure = [res_pressure; {string(M), height(Tk), t, df, p}];
end


T_pressure = T;
print_effectsizes(T_pressure, 'MVPA (Pressure)');

% Plot MVPAs

get_mvpaplots(lo)


%% ROI Supplementary: Intensity Encoding (Tables S2.1, S2.2)

% Paired t-tests: Low vs. High intensity within each ROI and modality,
% separately for CBP and Controls

dropmvpa = {'general','sound','mechanical','FM_PAIN','FM_MSS'};
nb_roi = neural_baseline_clean(~ismember(neural_baseline_clean.measure, dropmvpa), :);
allROIs_supp = cellstr(unique(nb_roi.measure));

for grp = {'CBP','HC'}
    grpStr = grp{1};
    fprintf('\n==== INTENSITY ENCODING: %s (Tables S2.1/S2.2) ====\n', grpStr);
    fprintf('%-16s | %-10s | %6s %6s | %8s %8s | %8s %8s %10s %6s\n', ...
        'ROI','Modality','nLow','nHigh','M_Low','M_High','g','p','q(BH)','');
    fprintf('%s\n', repmat('-',1,100));

    % Collect p-values for FDR across all ROI×Modality within this group
    pvals = [];
    testInfo = {};

    for mod = {'Sound','Pressure'}
        modStr = mod{1};
        for ri = 1:numel(allROIs_supp)
            R = allROIs_supp{ri};

            Tk = nb_roi(nb_roi.measure==R & nb_roi.modality==modStr & nb_roi.group==grpStr, :);
            if isempty(Tk), continue; end

            % Get paired data: match subjects across Low and High
            T_lo = Tk(Tk.intensity=="Low", {'subID','value'});
            T_hi = Tk(Tk.intensity=="High", {'subID','value'});
            [~, ia, ib] = intersect(T_lo.subID, T_hi.subID);

            xLo = T_lo.value(ia);
            xHi = T_hi.value(ib);

            % Remove pairs with NaN
            valid = ~isnan(xLo) & ~isnan(xHi);
            xLo = xLo(valid);  xHi = xHi(valid);

            if numel(xLo) < 3
                pvals(end+1) = NaN;
                testInfo{end+1} = {R, modStr, numel(xLo), numel(xHi), NaN, NaN, NaN, NaN};
                continue;
            end

            [~,p_tt] = ttest(xLo, xHi);
            ef = mes(xLo, xHi, 'hedgesg');

            pvals(end+1) = p_tt;
            testInfo{end+1} = {R, modStr, numel(xLo), numel(xHi), ...
                mean(xLo,'omitnan'), mean(xHi,'omitnan'), ef.hedgesg, p_tt};
        end
    end

    % FDR correction (Benjamini-Hochberg) across all tests within this group
    qvals = nan(size(pvals));
    keep = ~isnan(pvals);
    if any(keep), [~,~,~,qvals(keep)] = fdr_bh(pvals(keep)); end

    % Print
    for ti = 1:numel(testInfo)
        info = testInfo{ti};
        if isnan(info{7})
            fprintf('%-16s | %-10s | %6d %6d | %8s %8s | %8s %8s %10s %6s\n', ...
                info{1}, info{2}, info{3}, info{4}, 'NA','NA','NA','NA','NA','');
        else
            fprintf('%-16s | %-10s | %6d %6d | %8.2f %8.2f | %+8.2f %8.3g %10.3g %6s\n', ...
                info{1}, info{2}, info{3}, info{4}, info{5}, info{6}, ...
                info{7}, info{8}, qvals(ti), pstars(qvals(ti)));
        end
    end
end


%% ROI Supplementary: Modality Specificity (Tables S2.3, S2.4)

% Paired t-tests: Sound vs. Pressure within each ROI and intensity level,
% separately for CBP and Controls

for grp = {'CBP','HC'}
    grpStr = grp{1};
    fprintf('\n==== MODALITY SPECIFICITY: %s (Tables S2.3/S2.4) ====\n', grpStr);
    fprintf('%-16s | %-10s | %6s %6s | %8s %8s | %8s %8s %10s %6s\n', ...
        'ROI','Intensity','nSnd','nPrs','M_Sound','M_Press','g','p','q(BH)','');
    fprintf('%s\n', repmat('-',1,100));

    % Collect p-values for FDR across all ROI×Intensity within this group
    pvals = [];
    testInfo = {};

    for intLev = {'Low','High'}
        intStr = intLev{1};
        for ri = 1:numel(allROIs_supp)
            R = allROIs_supp{ri};

            Tk = nb_roi(nb_roi.measure==R & nb_roi.intensity==intStr & nb_roi.group==grpStr, :);
            if isempty(Tk), continue; end

            % Get paired data: match subjects across Sound and Pressure
            T_snd = Tk(Tk.modality=="Sound", {'subID','value'});
            T_prs = Tk(Tk.modality=="Pressure", {'subID','value'});
            [~, ia, ib] = intersect(T_snd.subID, T_prs.subID);

            xSnd = T_snd.value(ia);
            xPrs = T_prs.value(ib);

            % Remove pairs with NaN
            valid = ~isnan(xSnd) & ~isnan(xPrs);
            xSnd = xSnd(valid);  xPrs = xPrs(valid);

            if numel(xSnd) < 3
                pvals(end+1) = NaN;
                testInfo{end+1} = {R, intStr, numel(xSnd), numel(xPrs), NaN, NaN, NaN, NaN};
                continue;
            end

            [~,p_tt] = ttest(xSnd, xPrs);
            ef = mes(xSnd, xPrs, 'hedgesg');

            pvals(end+1) = p_tt;
            testInfo{end+1} = {R, intStr, numel(xSnd), numel(xPrs), ...
                mean(xSnd,'omitnan'), mean(xPrs,'omitnan'), ef.hedgesg, p_tt};
        end
    end

    % FDR correction (Benjamini-Hochberg) across all tests within this group
    qvals = nan(size(pvals));
    keep = ~isnan(pvals);
    if any(keep), [~,~,~,qvals(keep)] = fdr_bh(pvals(keep)); end

    % Print
    for ti = 1:numel(testInfo)
        info = testInfo{ti};
        if isnan(info{7})
            fprintf('%-16s | %-10s | %6d %6d | %8s %8s | %8s %8s %10s %6s\n', ...
                info{1}, info{2}, info{3}, info{4}, 'NA','NA','NA','NA','NA','');
        else
            fprintf('%-16s | %-10s | %6d %6d | %8.2f %8.2f | %+8.2f %8.3g %10.3g %6s\n', ...
                info{1}, info{2}, info{3}, info{4}, info{5}, info{6}, ...
                info{7}, info{8}, qvals(ti), pstars(qvals(ti)));
        end
    end
end


%% Neural Baseline Plot Statistics (with significance stars)

T = neural_baseline_clean;

% Ensure text columns are string/categorical for robust handling
for vn = ["group","measure","modality","intensity"]
    if ~iscategorical(T.(vn))
        T.(vn) = categorical(string(T.(vn)));
    end
end

% Explicitly set intensity order: Low before High
if ismember('Low', categories(T.intensity)) && ismember('High', categories(T.intensity))
    T.intensity = reordercats(T.intensity, {'Low','High'});
end

% Define combinations
allCombos = unique(T(:, {'measure','modality','intensity'}));

% Prepare results table
Results = table( ...
    strings(0,1), strings(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), ...
    NaN(0,1), NaN(0,1), NaN(0,1), strings(0,1), ...
    'VariableNames', { ...
        'measure','modality','intensity', ...
        'n_CBP','mean_CBP','sd_CBP', ...
        'n_HC','mean_HC','sd_HC', ...
        'tstat','df','pval','sigStars'});

for i = 1:height(allCombos)
    m  = allCombos.measure(i);
    mo = allCombos.modality(i);
    in = allCombos.intensity(i);

    % Slice data for this combo
    rows = T.measure == m & T.modality == mo & T.intensity == in;

    x = T.value(rows & T.group == 'CBP');
    y = T.value(rows & T.group == 'HC');

    % Skip if one group missing or too few obs
    if numel(x) < 2 || numel(y) < 2
        continue
    end

    % Descriptive stats
    nC = numel(x);  nH = numel(y);
    mC = mean(x);   mH = mean(y);
    sC = std(x, 0); sH = std(y, 0);

    % Welch’s t-test (unequal variances)
    [~, p, ~, stats] = ttest2(x, y, 'Vartype', 'unequal');

    % Round
    mC = round(mC, 2);
    mH = round(mH, 2);
    sC = round(sC, 2);
    sH = round(sH, 2);
    tval = round(stats.tstat, 2);
    df   = round(stats.df, 2);
    p    = round(p, 3);

    % --- Significance stars ---
    if p < 0.001
        stars = "***";
    elseif p < 0.01
        stars = "**";
    elseif p < 0.05
        stars = "*";
    else
        stars = "";
    end

    % Append to results
    Results = [Results; {
        string(m), string(mo), string(in), ...
        nC, mC, sC, ...
        nH, mH, sH, ...
        tval, df, p, string(stars)}];
end

% Sort with "Low" before "High"
Results.intensity = categorical(Results.intensity, {'Low','High'}, 'Ordinal', true);
Results = sortrows(Results, {'measure','modality','intensity'});

% Display
disp(Results)

% Optional: save to file
% writetable(Results, 'ttest_results_with_stars.csv');

%% Brain - Behavioural Correlations

fprintf('\n==== BRAIN-BEHAVIOR CORRELATIONS ====\n');

% Subset metadata for matching
M = d.metadata;
M = M(M.time == 1, :);                           % baseline only
isCBP = ismember(M.group, [1 2 3]);
M_cbp = M(isCBP, :);

% -- A1 vs. auditory unpleasantness (separate by group and intensity) --
fprintf('\n--- A1 vs. Auditory Unpleasantness ---\n');
fprintf('%-8s %-6s %8s %10s %6s\n', 'Group','Int','r','p','');
fprintf('%s\n', repmat('-',1,40));

for grp = {'HC','CBP'}
    grpStr = grp{1};
    if strcmp(grpStr,'HC')
        grpMask = ~isCBP;
    else
        grpMask = isCBP;
    end
    Mg = M(grpMask, :);

    % merge with A1 neural values
    for intPair = {{'Low','acute_mean_sound_lo','all_s_l'}, ...
                   {'High','acute_mean_sound_hi','all_s_h'}}
        ip = intPair{1};
        intLabel = ip{1}; behVar = ip{2}; neuField = ip{3};

        % Get A1 values for this group at baseline
        if strcmp(grpStr,'HC')
            a1_vals = lo.roi.HC.S1.A1.(neuField);
            ids_neural = d.HC.S1.id;
        else
            a1_vals = [lo.roi.G1.S1.A1.(neuField); lo.roi.G2.S1.A1.(neuField); lo.roi.G3.S1.A1.(neuField)];
            ids_neural = [d.G1.S1.id; d.G2.S1.id; d.G3.S1.id];
        end

        % Match by ID
        [~, ia, ib] = intersect(Mg.id, ids_neural);
        beh = Mg.(behVar)(ia);
        neu = a1_vals(ib);
        valid = ~isnan(beh) & ~isnan(neu);
        [r, p] = corr(beh(valid), neu(valid));
        fprintf('%-8s %-6s %8.2f %10.3g %6s\n', grpStr, intLabel, r, p, pstars(p));
    end
end

% -- mPFC vs. spontaneous pain (CBP only, low sound) --
fprintf('\n--- mPFC (Low Sound) vs. Spontaneous Pain (CBP only) ---\n');

mpfc_vals = [lo.roi.G1.S1.mPFC.all_s_l; lo.roi.G2.S1.mPFC.all_s_l; lo.roi.G3.S1.mPFC.all_s_l];
ids_neural = [d.G1.S1.id; d.G2.S1.id; d.G3.S1.id];

[~, ia, ib] = intersect(M_cbp.id, ids_neural);
beh = M_cbp.spon_pain_rating_mean(ia);
neu = mpfc_vals(ib);
valid = ~isnan(beh) & ~isnan(neu);
[r, p] = corr(beh(valid), neu(valid));
fprintf('  r(%d) = %.2f, p = %.3g%s\n', sum(valid)-2, r, p, pstars(p));

% -- Dorsal Insula vs. last-week back pain (CBP only, low sound) --
fprintf('\n--- Dorsal Insula (Low Sound) vs. Last-Week Pain (CBP only) ---\n');

dins_vals = [lo.roi.G1.S1.m_dorsal_insula.all_s_l; lo.roi.G2.S1.m_dorsal_insula.all_s_l; lo.roi.G3.S1.m_dorsal_insula.all_s_l];

[~, ia, ib] = intersect(M_cbp.id, ids_neural);
beh = M_cbp.pain_avg(ia);
neu = dins_vals(ib);
valid = ~isnan(beh) & ~isnan(neu);
[r, p] = corr(beh(valid), neu(valid));
fprintf('  r(%d) = %.2f, p = %.3g%s\n', sum(valid)-2, r, p, pstars(p));


%% Classification

% LASSO-regularized logistic regression with PCA (95% variance) and SMOTE
% 10-fold stratified cross-validation, 1000 iterations
% Three classifiers: behavioral, neural, combined
% (Reproduces run_classification.m / run_pca_classification.m)

nIter  = 1000;
nFolds = 10;

% ---- Prepare feature matrices using neural_baseline table ----

% Neural features: low-sound ROI averages + MVPA pattern expressions
% 15 features: 4 MVPA (FM_PAIN, FM_MSS, general, sound) + 11 ROIs
% NO mechanical pattern (matches get_classification_data 'sound_lo')
% Use neural_baseline (without outlier removal) for consistency with original

nb = neural_baseline(neural_baseline.modality=="Sound" & ...
     neural_baseline.intensity=="Low" & ...
     neural_baseline.measure ~= "mechanical", :);   % exclude mechanical
nb_wide = unstack(nb(:,{'subID','measure','value'}), 'value', 'measure');
nb_wide = sortrows(nb_wide, 'subID');

% Behavioral features: unpleasantness ratings for all 4 conditions
M_class = d.metadata;
M_class = M_class(M_class.time == 1, :);        % baseline only
M_class = sortrows(M_class, 'id');

% Align neural and behavioral by subject ID
[~, ia, ib] = intersect(M_class.id, nb_wide.subID);
M_class = M_class(ia, :);
nb_wide = nb_wide(ib, :);

X_beh  = [M_class.acute_mean_thumb_lo, M_class.acute_mean_thumb_hi, ...
          M_class.acute_mean_sound_lo, M_class.acute_mean_sound_hi];
X_neur = table2array(nb_wide(:, 2:end));         % drop subID column
X_comb = [X_beh, X_neur];

% Labels: CBP = 1, HC = 0
isCBP   = ismember(M_class.group, [1 2 3]);
y_class = double(isCBP);

fprintf('\nClassification features: %d behavioral, %d neural, %d combined\n', ...
    size(X_beh,2), size(X_neur,2), size(X_comb,2));
fprintf('Subjects: %d CBP, %d HC\n', sum(y_class==1), sum(y_class==0));

% ---- Run classifiers ----

classifiers = {'Behavioral','Neural','Combined'};
X_all = {X_beh, X_neur, X_comb};

fprintf('\n==== CLASSIFICATION: CBP vs. CONTROLS ====\n');
fprintf('%-12s | %8s ± %5s | %8s ± %5s | %8s ± %5s\n', ...
    'Classifier','AUC','sd','Sens','sd','Spec','sd');
fprintf('%s\n', repmat('-',1,72));

for ci = 1:numel(classifiers)
    x = X_all{ci};
    y = y_class;

    meanAUC  = zeros(nIter, 1);
    meanSens = zeros(nIter, 1);
    meanSpec = zeros(nIter, 1);

    for iteration = 1:nIter
        % PCA on full dataset (before CV split, matching original)
        [coeff, score, ~, ~, explained] = pca(x);

        % Select components that explain 95% of the variance
        explainedVariance = 0;
        numComponents = 0;
        while explainedVariance < 95
            numComponents = numComponents + 1;
            explainedVariance = sum(explained(1:numComponents));
        end
        x_pca = score(:, 1:numComponents);

        % Stratified 10-fold cross-validation
        cv = cvpartition(y, 'KFold', nFolds, 'Stratify', true);
        sensitivity = zeros(nFolds, 1);
        specificity = zeros(nFolds, 1);
        AUC = zeros(nFolds, 1);

        for i = 1:nFolds
            trainIdx = training(cv, i);
            X_train = x_pca(trainIdx, :);
            y_train = y(trainIdx);
            testIdx = test(cv, i);
            X_test = x_pca(testIdx, :);
            y_test = y(testIdx);

            % Handle class imbalance using SMOTE
            [X_train_balanced, y_train_balanced, ~, ~] = smote(X_train, [], 'Class', y_train);

            % Fit regularized logistic regression using lassoglm
            [B, FitInfo] = lassoglm(X_train_balanced, y_train_balanced, 'binomial', 'CV', 10);
            idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
            B0 = FitInfo.Intercept(idxLambdaMinDeviance);
            coef = [B0; B(:, idxLambdaMinDeviance)];

            % Predict class probabilities for test set
            X_test_intercept = [ones(size(X_test,1),1) X_test];
            scores = glmval(coef, X_test_intercept, 'logit', 'constant', 'off');

            % Calculate performance metrics using optimal ROC threshold
            [Xroc, Yroc, ~, AUC(i)] = perfcurve(y_test, scores, 1);
            [~, idx] = max(Yroc - Xroc);
            sensitivity(i) = Yroc(idx);
            specificity(i) = 1 - Xroc(idx);
        end

        % Store mean performance metrics for this iteration
        meanAUC(iteration)  = mean(AUC);
        meanSens(iteration) = mean(sensitivity);
        meanSpec(iteration) = mean(specificity);
    end

    fprintf('%-12s | %8.3f ± %.3f | %8.3f ± %.3f | %8.3f ± %.3f\n', ...
        classifiers{ci}, mean(meanAUC), std(meanAUC), ...
        mean(meanSens), std(meanSens), mean(meanSpec), std(meanSpec));
end


%% Exploratory Whole Brain Grey Matter Voxelwise Analysis

conds = {'all_s_l','all_s_h','all_t_l','all_t_h'};   % Sound Low/High, Pressure Low/High
alpha  = 0.001;        % voxelwise p for display (uncorrected)
kmin   = 10;           % min cluster size for display

for c = 1:numel(conds)
    cond = conds{c};
    fprintf('\n=== Condition: %s ===\n', cond);

    % -- 1) Gather fmri_data for CBP (G1+G2+G3) and HC at baseline (S1)
    Dcbp = [lo.G1.S1.(cond), lo.G2.S1.(cond), lo.G3.S1.(cond)];  % concatenate images (subjects)
    Dhc  =  lo.HC.S1.(cond);

    % -- 2) Concatenate CBP and HC into one dataset
    dat = [Dcbp, Dhc];   % fmri_data horzcat (images are columns)

    % Counts
    nCBP = size(Dcbp.dat, 2);
    nHC  = size(Dhc.dat,  2);
    nTot = nCBP + nHC;

    % -- 3) Smooth
    dat_smooth = preprocess(dat, 'smooth', 6);

    % -- 4) Weighted group regressor (sums to 0 so it's a clean contrast)
    g = [ ones(nCBP,1) * ( 1/nCBP );    ...  % CBP
         -ones(nHC,1)  * ( 1/nHC  ) ];      % HC
    X = [g, ones(nTot,1)];                  % [GroupContrast, Intercept]
    varnames = {'Group_CBP_vs_HC','Intercept'};

    % Quick check (should be ~0)
    fprintf('  mean(group regressor) = %.3e\n', mean(g));

    % -- 5) Mask to gray matter
    gm = fmri_data(which('gray_matter_mask.nii'));
    gm.dat = gm.dat > 0.3;
    dat_mask = apply_mask(dat_smooth, gm);

    % -- 6) Run regression (robust, uncorrected display thresholds)
    dat_mask.X = X;
    out = regress(dat_mask, 'variable_names', varnames, alpha, 'unc', 'k', kmin, 'robust');

    % -- 7) Pull the Group effect (1st regressor) and show
    Tmap = select_one_image(out.t, 1);     % 1 = Group_CBP_vs_HC
    figure('Color','w'); orthviews(Tmap);  title(sprintf('%s : CBP > HC (weighted)', cond));

    % Optional: slices montage
    figure('Color','w');
    montage(Tmap, 'trans', 'full');
    sgtitle(sprintf('CBP vs HC (weighted) — %s', cond));

    % -- 8) Table of clusters with (rough) labels
    reg = region(out.t);    % uses current threshold set by regress
    tbl = table(reg);
    disp(tbl);

    % ---- Optional saves ----
    % write(Tmap, 'fname', sprintf('tmap_CBP_vs_HC_%s.nii', cond), 'thresh');  % thresholded
    % write(out.t, 'fname', sprintf('tall_CBP_vs_HC_%s.nii', cond));           % full t-stack
end



%% Neural Longitudinal: LMM SOUND

% Fit & test treatment effects with LMM (random slopes & intercepts)

[results, models] = analyze_neural_longitudinal_sound(neural_longitudinal_sound);

%% Treatmenteffects
 
% Behavioral
plot_auditory_treatmenteffects_collapsed(d);

% Neural
T = neural_longitudinal_sound;
f = plot_mpfc_treatmenteffects_sound(T);



%% Predict treatment changes

% Behavioral
predict_treatment_changes(behavioral_longitudinal,d)

% Neural
predict_treatment_changes_neural(neural_longitudinal_sound, d)

% OUtcome: 