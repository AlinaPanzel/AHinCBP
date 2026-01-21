function plot_predictive_treatment_relationship

%% LOW
% 
T = predictTbl;
mdl = fitlm(T, 'pain_post ~ pain_pre + MSS_low_pre + MSS_low_pre:isPRT');
disp(mdl)
figure('Color','w'); hold on;

% Generate MSS range
x = linspace(min(T.MSS_low_pre), max(T.MSS_low_pre), 100)';

T.GroupCount = [];
% Hold pain_pre constant at mean
pain_pre_mean = mean(T.pain_pre,'omitnan');

% --- Raw data (faint, like ggplot background) ---
idx0 = T.isPRT == 0;
idx1 = T.isPRT == 1;

scatter(T.MSS_low_pre(idx0), T.pain_post(idx0), ...
    18, [0.4 0.6 1], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

scatter(T.MSS_low_pre(idx1), T.pain_post(idx1), ...
    18, [0.7 0.2 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

% --- Prediction lines + CI ---
for g = 0:1
    newTbl = table();
    newTbl.pain_pre = repmat(pain_pre_mean, numel(x), 1);
    newTbl.MSS_low_pre = x;
    newTbl.isPRT = repmat(g, numel(x), 1);

    [yhat, yCI] = predict(mdl, newTbl);

    if g == 0
        col = [0.2 0.5 0.9];   % control (blue)
        lbl = 'Control';
    else
        col = [0.5 0 0.6];     % PRT (purple)
        lbl = 'PRT';
    end

    % Confidence band
    fill([x; flipud(x)], ...
         [yCI(:,1); flipud(yCI(:,2))], ...
         col, 'FaceAlpha',0.15, 'EdgeColor','none');

    % Regression line
    plot(x, yhat, 'Color', col, 'LineWidth', 3, ...
         'DisplayName', lbl);
end

% --- Vertical zero-line if you want the same reference style ---
xline(0,'--k','LineWidth',1);

% --- Formatting to match your example ---
xlabel('Baseline Auditory Unpleasantness (Low Intensity)');
ylabel('Post-Treatment Pain');
title('Auditory Hyperresponsivity × Treatment Interaction');
grid on;
box on;
set(gca,'FontSize',12);
%legend('Location','northwest');

hold off;



%% HIGH

T = predictTbl;
mdl = fitlm(T, 'pain_post ~ pain_pre + MSS_high_pre + MSS_high_pre:isPRT');
disp(mdl)
figure('Color','w'); hold on;

% Generate MSS range
x = linspace(min(T.MSS_high_pre), max(T.MSS_high_pre), 100)';

% Hold pain_pre constant at mean
pain_pre_mean = mean(T.pain_pre,'omitnan');

% --- Raw data ---
idx0 = T.isPRT == 0;
idx1 = T.isPRT == 1;

scatter(T.MSS_high_pre(idx0), T.pain_post(idx0), ...
    18, [0.4 0.6 1], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

scatter(T.MSS_high_pre(idx1), T.pain_post(idx1), ...
    18, [0.7 0.2 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

% --- Prediction lines + CI ---
for g = 0:1
    newTbl = table();
    newTbl.pain_pre = repmat(pain_pre_mean, numel(x), 1);
    newTbl.MSS_high_pre = x;
    newTbl.isPRT = repmat(g, numel(x), 1);

    [yhat, yCI] = predict(mdl, newTbl);

    if g == 0
        col = [0.2 0.5 0.9];   % control (blue)
        lbl = 'Control';
    else
        col = [0.5 0 0.6];     % PRT (purple)
        lbl = 'PRT';
    end

    % Confidence band
    fill([x; flipud(x)], ...
         [yCI(:,1); flipud(yCI(:,2))], ...
         col, 'FaceAlpha',0.15, 'EdgeColor','none');

    % Regression line
    plot(x, yhat, 'Color', col, 'LineWidth', 3, ...
         'DisplayName', lbl);
end

% --- Vertical zero-line if you want the same reference style ---
xline(0,'--k','LineWidth',1);

% --- Formatting to match your example ---
xlabel('Baseline Auditory Unpleasantness (High Intensity)');
ylabel('Post-Treatment Pain');
title('Auditory Hyperresponsivity × Treatment Interaction');
grid on;
box on;
set(gca,'FontSize',12);
%legend('Location','northwest');

hold off;

end
