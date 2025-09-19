%% Auditory Hypersensitivity in Chronic Back Pain: A Randomized Controlled Trial of Pain Reprocessing Therapy


% Author: Alina Panzel
% Last Date of Changes: 19.09.2025

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



% PLOT AUDITORY 

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Auditory");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [-2 7];

% Plot data for different regions
plot_allgoodplot(roi.A1, [1 2 3 4], 'Auditory Cortex','roi', yLimit, 3, true);
plot_allgoodplot(roi.MGN, [5 6 7 8], 'Medial Geniculate Nucleus','roi', yLimit, 3, true);
plot_allgoodplot(roi.IC, [9 10 11 12], 'Inferior Colliculus','roi', yLimit, 3, true);






% -------- Behavioural Longitudinal Analysis ----------

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


%% Neural Baseline: LMM ROI

fprintf('\n==== LMM RESULTS SOUND ====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

res_sound = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','df1','df2','p_group'});

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


    res_sound = [res_sound; {string(R), height(Tk), T, DF1, DF2, p}];
end


fprintf('\n==== EFFECT SIZES SOUND (ROIs) ====\n');
fprintf('%-16s | %s | %s\n','ROI','Low Intensity','High Intensity');
fprintf('%-16s | %s | %s\n',' ', ...
    'HC (M±SD), CBP (M±SD), g [95%% CI]', ...
    'HC (M±SD), CBP (M±SD), g [95%% CI]');
fprintf('%s\n', repmat('-',1,120));

for i = 1:numel(allROIs)
    R  = allROIs{i};
    Tk = neural_baseline_roi_sound(neural_baseline_roi_sound.measure==R, :);
    if isempty(Tk), continue; end

    % --- Low intensity ---
    xL = Tk.value(Tk.group=="HC"  & Tk.intensity=="Low");
    yL = Tk.value(Tk.group=="CBP" & Tk.intensity=="Low");
    mHCL = mean(xL,'omitnan'); sHCL = std(xL,'omitnan');
    mCBPL = mean(yL,'omitnan'); sCBPL = std(yL,'omitnan');
    efL = mes(xL,yL,'hedgesg');                 % Hedges' g (MES toolbox)
    gL  = efL.hedgesg;  ciL = efL.hedgesgCi(:)';

    % --- High intensity ---
    xH = Tk.value(Tk.group=="HC"  & Tk.intensity=="High");
    yH = Tk.value(Tk.group=="CBP" & Tk.intensity=="High");
    mHCH = mean(xH,'omitnan'); sHCH = std(xH,'omitnan');
    mCBPH = mean(yH,'omitnan'); sCBPH = std(yH,'omitnan');
    efH = mes(xH,yH,'hedgesg');
    gH  = efH.hedgesg;  ciH = efH.hedgesgCi(:)';

    fprintf('%-16s | HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f] | HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f]\n', ...
        R, mHCL,sHCL,mCBPL,sCBPL,gL,ciL(1),ciL(2), ...
           mHCH,sHCH,mCBPH,sCBPH,gH,ciH(1),ciH(2));
end


fprintf('\n==== LMM RESULTS PRESSURE====\n');
fprintf('%-16s | %4s %4s | %8s %8s | %8s %8s | %10s\n', ...
    'ROI','nHC','nCBP','HC_Low','HC_High','CBP_Low','CBP_High','Group T [DF] (p)');
fprintf('%s\n', repmat('-',1,90));

dropmvpa = {'general','sound','mechanical','FM_PAIN','FM_MSS'};

res_pressure = table(string.empty, zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','N','T_group','df1','df2','p_group'});

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
    res_pressure = [res_pressure; {string(R), height(Tk), T, DF1, DF2, p}];
end


fprintf('\n==== EFFECT SIZES PRESSURE (ROIs) ====\n');
fprintf('%-16s | %s | %s\n','ROI','Low Intensity','High Intensity');
fprintf('%-16s | %s | %s\n',' ', ...
    'HC (M±SD), CBP (M±SD), g [95%% CI]', ...
    'HC (M±SD), CBP (M±SD), g [95%% CI]');
fprintf('%s\n', repmat('-',1,120));

for i = 1:numel(allROIs)
    R  = allROIs{i};
    Tk = neural_baseline_roi_pressure(neural_baseline_roi_pressure.measure==R, :);
    if isempty(Tk), continue; end

    % --- Low intensity ---
    xL = Tk.value(Tk.group=="HC"  & Tk.intensity=="Low");
    yL = Tk.value(Tk.group=="CBP" & Tk.intensity=="Low");
    mHCL = mean(xL,'omitnan'); sHCL = std(xL,'omitnan');
    mCBPL = mean(yL,'omitnan'); sCBPL = std(yL,'omitnan');
    efL = mes(xL,yL,'hedgesg');                 % Hedges' g (MES toolbox)
    gL  = efL.hedgesg;  ciL = efL.hedgesgCi(:)';

    % --- High intensity ---
    xH = Tk.value(Tk.group=="HC"  & Tk.intensity=="High");
    yH = Tk.value(Tk.group=="CBP" & Tk.intensity=="High");
    mHCH = mean(xH,'omitnan'); sHCH = std(xH,'omitnan');
    mCBPH = mean(yH,'omitnan'); sCBPH = std(yH,'omitnan');
    efH = mes(xH,yH,'hedgesg');
    gH  = efH.hedgesg;  ciH = efH.hedgesgCi(:)';

    fprintf('%-16s | HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f] | HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f]\n', ...
        R, mHCL,sHCL,mCBPL,sCBPL,gL,ciL(1),ciL(2), ...
           mHCH,sHCH,mCBPH,sCBPH,gH,ciH(1),ciH(2));
end

%% Neural Baseline: LMM MVPA
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


fprintf('\n==== EFFECT SIZES SOUND (MVPAs) ====\n');
fprintf('%-16s | %10s | %10s\n','MVPA','g_Low [95% CI]','g_High [95% CI]');
fprintf('%s\n', repmat('-',1,46));

for i = 1:numel(allMVPAs)
    M  = allMVPAs{i};
    Tk = T(T.measure==M,:);
    if isempty(Tk), continue; end
 % Low intensity
    xL = Tk.value(Tk.group=="HC"  & Tk.intensity=="Low");
    yL = Tk.value(Tk.group=="CBP" & Tk.intensity=="Low");
    efL = mes(xL,yL,'hedgesg');
    gL  = efL.hedgesg;
    ciL = efL.hedgesgCi(:)';
    mHC_L = mean(xL,'omitnan'); sHC_L = std(xL,'omitnan');
    mCBP_L = mean(yL,'omitnan'); sCBP_L = std(yL,'omitnan');

    % High intensity
    xH = Tk.value(Tk.group=="HC"  & Tk.intensity=="High");
    yH = Tk.value(Tk.group=="CBP" & Tk.intensity=="High");
    efH = mes(xH,yH,'hedgesg');
    gH  = efH.hedgesg;
    ciH = efH.hedgesgCi(:)';
    mHC_H = mean(xH,'omitnan'); sHC_H = std(xH,'omitnan');
    mCBP_H = mean(yH,'omitnan'); sCBP_H = std(yH,'omitnan');

    % Print one row: effect sizes + descriptive stats
    fprintf('%-16s | HC %.2f±%.2f, CBP %.2f±%.2f, g=%.2f [%+.2f, %+.2f] | HC %.2f±%.2f, CBP %.2f±%.2f, g=%.2f [%+.2f, %+.2f]\n', ...
        M, mHC_L, sHC_L, mCBP_L, sCBP_L, gL, ciL(1), ciL(2), ...
           mHC_H, sHC_H, mCBP_H, sCBP_H, gH, ciH(1), ciH(2));
end



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

fprintf('\n==== EFFECT SIZES PRESSURE (MVPAs) ====\n');
fprintf('%-16s | %24s | %24s\n','MVPA','Low Intensity','High Intensity');
fprintf('%-16s | %9s %9s | %9s %9s\n',' ', 'HC (M±SD)','CBP (M±SD)','HC (M±SD)','CBP (M±SD)');
fprintf('%s\n', repmat('-',1,80));

for i = 1:numel(allMVPAs)
    M  = allMVPAs{i};
    Tk = T(T.measure==M,:);
    if isempty(Tk), continue; end

    % Low intensity
    xL = Tk.value(Tk.group=="HC"  & Tk.intensity=="Low");
    yL = Tk.value(Tk.group=="CBP" & Tk.intensity=="Low");
    efL = mes(xL,yL,'hedgesg');
    gL  = efL.hedgesg;
    ciL = efL.hedgesgCi(:)';
    mHC_L = mean(xL,'omitnan'); sHC_L = std(xL,'omitnan');
    mCBP_L = mean(yL,'omitnan'); sCBP_L = std(yL,'omitnan');

    % High intensity
    xH = Tk.value(Tk.group=="HC"  & Tk.intensity=="High");
    yH = Tk.value(Tk.group=="CBP" & Tk.intensity=="High");
    efH = mes(xH,yH,'hedgesg');
    gH  = efH.hedgesg;
    ciH = efH.hedgesgCi(:)';
    mHC_H = mean(xH,'omitnan'); sHC_H = std(xH,'omitnan');
    mCBP_H = mean(yH,'omitnan'); sCBP_H = std(yH,'omitnan');

    % Print one row: effect sizes + descriptive stats
    fprintf('%-16s | HC %.2f±%.2f, CBP %.2f±%.2f, g=%.2f [%+.2f, %+.2f] | HC %.2f±%.2f, CBP %.2f±%.2f, g=%.2f [%+.2f, %+.2f]\n', ...
        M, mHC_L, sHC_L, mCBP_L, sCBP_L, gL, ciL(1), ciL(2), ...
           mHC_H, sHC_H, mCBP_H, sCBP_H, gH, ciH(1), ciH(2));
end


%% Classification

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
    models(end+1).measure = string(m); 
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




