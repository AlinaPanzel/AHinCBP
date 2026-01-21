%% mPFC Moderation Analysis - Fixed version

% 1. Extract mPFC baseline data (using contains to handle whitespace)
mPFC_baseline = neural_longitudinal_sound(contains(string(neural_longitudinal_sound.measure), 'mPFC') & ...
                                          neural_longitudinal_sound.timepoint == 1, :);

% 2. Pivot to wide format (separate columns for low and high intensity)
mPFC_low = mPFC_baseline(mPFC_baseline.intensity == 1, {'subID', 'group', 'value'});
mPFC_low.Properties.VariableNames{'value'} = 'mPFC_low_pre';

mPFC_high = mPFC_baseline(mPFC_baseline.intensity == 2, {'subID', 'value'});
mPFC_high.Properties.VariableNames{'value'} = 'mPFC_high_pre';

% Merge low and high
mPFC_wide = outerjoin(mPFC_low, mPFC_high, 'Keys', 'subID', 'MergeKeys', true);

fprintf('mPFC wide format: %d subjects\n', height(mPFC_wide));

% 3. Merge with behavioral data
% predictTbl has SubID (capital S), mPFC_wide has subID (lowercase)
mPFC_wide.Properties.VariableNames{'subID'} = 'SubID';

mergedTbl = innerjoin(predictTbl, mPFC_wide(:, {'SubID', 'mPFC_low_pre', 'mPFC_high_pre'}), 'Keys', 'SubID');
fprintf('Merged table: %d subjects\n', height(mergedTbl));

% 4. Correlate MSS with mPFC
[r_low, p_low] = corr(mergedTbl.MSS_low_pre, mergedTbl.mPFC_low_pre, 'rows', 'complete');
[r_high, p_high] = corr(mergedTbl.MSS_high_pre, mergedTbl.mPFC_high_pre, 'rows', 'complete');

fprintf('\n=== MSS vs mPFC Correlations ===\n');
fprintf('Low intensity:  r = %.3f, p = %.4f\n', r_low, p_low);
fprintf('High intensity: r = %.3f, p = %.4f\n', r_high, p_high);

% 5. Scatter plot
figure('Position', [100 100 900 400]);
subplot(1,2,1);
scatter(mergedTbl.MSS_low_pre, mergedTbl.mPFC_low_pre, 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('MSS Low (Auditory Unpleasantness)');
ylabel('mPFC Activity (Low Intensity)');
title(sprintf('Low Intensity: r = %.2f, p = %.3f', r_low, p_low));
lsline; grid on;

subplot(1,2,2);
scatter(mergedTbl.MSS_high_pre, mergedTbl.mPFC_high_pre, 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('MSS High (Auditory Unpleasantness)');
ylabel('mPFC Activity (High Intensity)');
title(sprintf('High Intensity: r = %.2f, p = %.3f', r_high, p_high));
lsline; grid on;

sgtitle('MSS (Behavioral) vs mPFC Activity (Neural)');
saveas(gcf, 'MSS_mPFC_correlation.png');
fprintf('Saved: MSS_mPFC_correlation.png\n');



%% PARALLEL MODERATION ANALYSIS: MSS vs mPFC as treatment moderators

% 1. Prepare merged data (already done above, but ensure we have it)
fprintf('=== PARALLEL MODERATION: MSS vs mPFC ===\n\n');
fprintf('Subjects in merged analysis: %d\n', height(mergedTbl));

% 2. Set up group coding (Placebo as reference)
group_cat = categorical(mergedTbl.group, [2 1 3], {'Placebo', 'PRT', 'UC'});
mergedTbl.group_cat = group_cat;

% 3. MSS-mPFC correlation (secondary - just to see relationship)
[r_low, p_low] = corr(mergedTbl.MSS_low_pre, mergedTbl.mPFC_low_pre, 'rows', 'complete');
fprintf('MSS-mPFC correlation (low): r = %.3f, p = %.4f\n\n', r_low, p_low);

%% MODEL 1: MSS as moderator (rerun for comparison)
fprintf('=== MODEL 1: MSS as Moderator ===\n');
mdl_MSS = fitlm(mergedTbl, 'pain_post ~ pain_pre + MSS_low_pre*group_cat');
disp(mdl_MSS.Coefficients);

%% MODEL 2: mPFC as moderator (NEW)
fprintf('\n=== MODEL 2: mPFC as Moderator ===\n');
mdl_mPFC = fitlm(mergedTbl, 'pain_post ~ pain_pre + mPFC_low_pre*group_cat');
disp(mdl_mPFC.Coefficients);

%% Compare model fits
fprintf('\n=== Model Comparison ===\n');
fprintf('MSS model R² = %.3f, Adjusted R² = %.3f\n', mdl_MSS.Rsquared.Ordinary, mdl_MSS.Rsquared.Adjusted);
fprintf('mPFC model R² = %.3f, Adjusted R² = %.3f\n', mdl_mPFC.Rsquared.Ordinary, mdl_mPFC.Rsquared.Adjusted);



%% MODERATION ANALYSIS: Insula subregions and A1

% Extract baseline neural measures (timepoint == 1)
T_neural = neural_longitudinal_sound;
baseline = T_neural(T_neural.timepoint == 1, :);

% Check available measures
disp('Available measures:');
disp(unique(baseline.measure));

%% Extract each ROI for low intensity (intensity == 1)
measures = {'A1', 'm_ventral_insula', 'm_dorsal_insula', 'm_posterior_insula'};
shortnames = {'A1', 'vIns', 'dIns', 'pIns'};

% Create table with all neural predictors
neural_wide = table();

for i = 1:length(measures)
    % Use contains to handle whitespace
    idx = contains(string(baseline.measure), measures{i}) & baseline.intensity == 1;
    tmp = baseline(idx, {'subID', 'value'});
    tmp.Properties.VariableNames{'value'} = [shortnames{i} '_low_pre'];
    
    if isempty(neural_wide)
        neural_wide = tmp;
    else
        neural_wide = outerjoin(neural_wide, tmp, 'Keys', 'subID', 'MergeKeys', true);
    end
end

% Add mPFC (already know it works)
idx = contains(string(baseline.measure), 'mPFC') & baseline.intensity == 1;
tmp = baseline(idx, {'subID', 'value'});
tmp.Properties.VariableNames{'value'} = 'mPFC_low_pre';
neural_wide = outerjoin(neural_wide, tmp, 'Keys', 'subID', 'MergeKeys', true);

fprintf('Neural predictors extracted for %d subjects\n', height(neural_wide));

%% Merge with behavioral data
neural_wide.Properties.VariableNames{'subID'} = 'SubID';
mergedTbl2 = innerjoin(predictTbl, neural_wide, 'Keys', 'SubID');

% Set up group coding (Placebo as reference)
mergedTbl2.group_cat = categorical(mergedTbl2.group, [2 1 3], {'Placebo', 'PRT', 'UC'});

fprintf('Merged: %d subjects for analysis\n\n', height(mergedTbl2));

%% Run moderation models for each neural ROI
fprintf('=================================================================\n');
fprintf('MODERATION ANALYSIS: Neural ROIs as Treatment Moderators\n');
fprintf('Reference group: Placebo\n');
fprintf('=================================================================\n\n');

roi_names = {'A1_low_pre', 'vIns_low_pre', 'dIns_low_pre', 'pIns_low_pre', 'mPFC_low_pre'};
roi_labels = {'A1 (Primary Auditory)', 'Ventral Insula', 'Dorsal Insula', 'Posterior Insula', 'mPFC'};

results = table();

for i = 1:length(roi_names)
    roi = roi_names{i};
    
    % Build formula
    formula = sprintf('pain_post ~ pain_pre + %s*group_cat', roi);
    
    % Fit model
    mdl = fitlm(mergedTbl2, formula);
    
    % Extract key statistics
    coef_tbl = mdl.Coefficients;
    
    % Find the PRT interaction term
    prt_int_idx = contains(coef_tbl.Properties.RowNames, 'group_cat_PRT') & ...
                  contains(coef_tbl.Properties.RowNames, roi);
    uc_int_idx = contains(coef_tbl.Properties.RowNames, 'group_cat_UC') & ...
                 contains(coef_tbl.Properties.RowNames, roi);
    
    % Store results
    results.ROI{i} = roi_labels{i};
    results.R2(i) = mdl.Rsquared.Ordinary;
    results.R2_adj(i) = mdl.Rsquared.Adjusted;
    
    % PRT interaction
    if any(prt_int_idx)
        results.PRT_beta(i) = coef_tbl.Estimate(prt_int_idx);
        results.PRT_p(i) = coef_tbl.pValue(prt_int_idx);
    end
    
    % UC interaction
    if any(uc_int_idx)
        results.UC_beta(i) = coef_tbl.Estimate(uc_int_idx);
        results.UC_p(i) = coef_tbl.pValue(uc_int_idx);
    end
    
    % Print full model
    fprintf('--- %s ---\n', roi_labels{i});
    disp(coef_tbl);
    fprintf('\n');
end

%% Summary table
fprintf('=================================================================\n');
fprintf('SUMMARY: Neural Moderation of Treatment Response\n');
fprintf('=================================================================\n');
results.PRT_sig = results.PRT_p < 0.05;
results.UC_sig = results.UC_p < 0.05;
disp(results);

%% Add MSS for comparison
fprintf('\n--- MSS (Behavioral - for comparison) ---\n');
mdl_MSS = fitlm(mergedTbl2, 'pain_post ~ pain_pre + MSS_low_pre*group_cat');
disp(mdl_MSS.Coefficients);
fprintf('R² = %.3f\n', mdl_MSS.Rsquared.Ordinary);