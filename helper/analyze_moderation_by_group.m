function analyze_moderation_by_group()
%% FULL MODEL - Three groups
% Group 1 = PRT, Group 2 = Placebo, Group 3 = Usual Care

T = predictTbl;

% Convert group to categorical with PRT (group 1) as reference
T.group = categorical(T.group, [1 2 3], {'PRT','Placebo','UC'});

%% LOW INTENSITY MODEL
mdl_low = fitlm(T, 'pain_post ~ pain_pre + MSS_low_pre + group + MSS_low_pre:group');
disp('=== LOW INTENSITY SOUNDS (PRT = reference) ===')
disp(mdl_low)

%% HIGH INTENSITY MODEL
mdl_high = fitlm(T, 'pain_post ~ pain_pre + MSS_high_pre + group + MSS_high_pre:group');
disp('=== HIGH INTENSITY SOUNDS (PRT = reference) ===')
disp(mdl_high)

%% FIGURE - LOW INTENSITY
figure('Color','w', 'Position', [100 100 600 500]); hold on;

% Colors
col_PRT = [0.5 0 0.6];       % purple
col_Placebo = [0.2 0.6 0.2]; % green
col_UC = [0.2 0.5 0.9];      % blue

% MSS range for predictions
x = linspace(min(T.MSS_low_pre), max(T.MSS_low_pre), 100)';
pain_pre_mean = mean(T.pain_pre, 'omitnan');

% --- Raw data points ---
idx_PRT = T.group == 'PRT';
idx_Placebo = T.group == 'Placebo';
idx_UC = T.group == 'UC';

scatter(T.MSS_low_pre(idx_UC), T.pain_post(idx_UC), 25, col_UC, 'filled', 'MarkerFaceAlpha', 0.3);
scatter(T.MSS_low_pre(idx_Placebo), T.pain_post(idx_Placebo), 25, col_Placebo, 'filled', 'MarkerFaceAlpha', 0.3);
scatter(T.MSS_low_pre(idx_PRT), T.pain_post(idx_PRT), 25, col_PRT, 'filled', 'MarkerFaceAlpha', 0.3);

% --- Prediction lines + CI for each group ---
groups = {'PRT', 'Placebo', 'UC'};
colors = {col_PRT, col_Placebo, col_UC};

h_lines = gobjects(3,1);

for g = 1:3
    newTbl = table();
    newTbl.pain_pre = repmat(pain_pre_mean, numel(x), 1);
    newTbl.MSS_low_pre = x;
    newTbl.group = repmat(categorical(groups(g), {'PRT','Placebo','UC'}), numel(x), 1);
    
    [yhat, yCI] = predict(mdl_low, newTbl);
    
    % Confidence band
    fill([x; flipud(x)], [yCI(:,1); flipud(yCI(:,2))], ...
         colors{g}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Regression line
    h_lines(g) = plot(x, yhat, 'Color', colors{g}, 'LineWidth', 2.5);
end

xlabel('Baseline Auditory Unpleasantness (Low Intensity)', 'FontSize', 12);
ylabel('Post-Treatment Pain (0-10)', 'FontSize', 12);
title('Treatment Moderation by Auditory Hyperresponsivity', 'FontSize', 14);
legend(h_lines, groups, 'Location', 'northeast');
grid on; box on;
set(gca, 'FontSize', 11);

hold off;

%% FIGURE - HIGH INTENSITY
figure('Color','w', 'Position', [100 100 600 500]); hold on;

x = linspace(min(T.MSS_high_pre), max(T.MSS_high_pre), 100)';

scatter(T.MSS_high_pre(idx_UC), T.pain_post(idx_UC), 25, col_UC, 'filled', 'MarkerFaceAlpha', 0.3);
scatter(T.MSS_high_pre(idx_Placebo), T.pain_post(idx_Placebo), 25, col_Placebo, 'filled', 'MarkerFaceAlpha', 0.3);
scatter(T.MSS_high_pre(idx_PRT), T.pain_post(idx_PRT), 25, col_PRT, 'filled', 'MarkerFaceAlpha', 0.3);

h_lines = gobjects(3,1);

for g = 1:3
    newTbl = table();
    newTbl.pain_pre = repmat(pain_pre_mean, numel(x), 1);
    newTbl.MSS_high_pre = x;
    newTbl.group = repmat(categorical(groups(g), {'PRT','Placebo','UC'}), numel(x), 1);
    
    [yhat, yCI] = predict(mdl_high, newTbl);
    
    fill([x; flipud(x)], [yCI(:,1); flipud(yCI(:,2))], ...
         colors{g}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    h_lines(g) = plot(x, yhat, 'Color', colors{g}, 'LineWidth', 2.5);
end

xlabel('Baseline Auditory Unpleasantness (High Intensity)', 'FontSize', 12);
ylabel('Post-Treatment Pain (0-10)', 'FontSize', 12);
title('Treatment Moderation by Auditory Hyperresponsivity', 'FontSize', 14);
legend(h_lines, groups, 'Location', 'northeast');
grid on; box on;
set(gca, 'FontSize', 11);

hold off;

end
