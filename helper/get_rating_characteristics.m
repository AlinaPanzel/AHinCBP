function baseline_ratings = get_rating_characteristics(d)

% Subset to CBP patients at baseline only
M = d.metadata;
M = M(M.time == 1, :);
isCBP = ismember(M.group, [1 2 3]);
M = M(isCBP, :);

ids   = M.id;
spon  = M.spon_pain_ratings;      % N x 8
unpl  = M.acute_pain_ratings;     % N x 10
thumb = M.acute_stim_is_thumb;    % N x 10 (0=sound, 1=thumb)
inten = M.acute_stim_intensity;   % N x 10 (1=low, 2=high)

% Spontaneous pain metrics
spon_mean = mean(spon, 2, 'omitnan');
spon_var  = var(spon, 0, 2, 'omitnan');

% Preallocate outputs
nSubjects = size(unpl,1);

sound_low_mean     = nan(nSubjects, 1);
sound_low_var      = nan(nSubjects, 1);
sound_high_mean    = nan(nSubjects, 1);
sound_high_var     = nan(nSubjects, 1);

press_low_mean     = nan(nSubjects, 1);
press_low_var      = nan(nSubjects, 1);
press_high_mean    = nan(nSubjects, 1);
press_high_var     = nan(nSubjects, 1);

% Loop per subject
for s = 1:nSubjects
    
    % Logical indices for each condition
    idx_sound_low   = (thumb(s,:) == 0) & (inten(s,:) == 1);
    idx_sound_high  = (thumb(s,:) == 0) & (inten(s,:) == 2);
    idx_press_low   = (thumb(s,:) == 1) & (inten(s,:) == 1);
    idx_press_high  = (thumb(s,:) == 1) & (inten(s,:) == 2);

    % Compute condition means / variances
    sound_low_mean(s)   = mean(unpl(s, idx_sound_low), 'omitnan');
    sound_low_var(s)    = var(unpl(s, idx_sound_low), 'omitnan');

    sound_high_mean(s)  = mean(unpl(s, idx_sound_high), 'omitnan');
    sound_high_var(s)   = var(unpl(s, idx_sound_high), 'omitnan');

    press_low_mean(s)   = mean(unpl(s, idx_press_low), 'omitnan');
    press_low_var(s)    = var(unpl(s, idx_press_low), 'omitnan');

    press_high_mean(s)  = mean(unpl(s, idx_press_high), 'omitnan');
    press_high_var(s)   = var(unpl(s, idx_press_high), 'omitnan');
end

% Build table
summaryTable = table(ids, ...
    spon_mean, spon_var, ...
    sound_low_mean, sound_low_var, ...
    sound_high_mean, sound_high_var, ...
    press_low_mean, press_low_var, ...
    press_high_mean, press_high_var);

disp(summaryTable);

baseline_ratings = summaryTable;

%% === CORRELATIONS (6 x 6) ===

% Extract variables
sp_mean   = summaryTable.spon_mean;
sp_var    = summaryTable.spon_var;
sound_low = summaryTable.sound_low_mean;
sound_high = summaryTable.sound_high_mean;
press_low  = summaryTable.press_low_mean;
press_high = summaryTable.press_high_mean;

% All variables in one matrix (columns = variables)
allVars = [sp_mean, sp_var, sound_low, sound_high, press_low, press_high];

% Names in same order as columns
allNames = {'Spont Mean', 'Spont Var', ...
            'Sound Low', 'Sound High', ...
            'Press Low', 'Press High'};

% Full 6x6 correlation and p-value matrices
[R, P] = corr(allVars, 'rows', 'pairwise');

%% === PRINT RESULTS (optional) ===
fprintf('\n=== FULL 6x6 CORRELATION MATRIX ===\n');
disp(array2table(R, 'VariableNames', allNames, 'RowNames', allNames));

% Round full correlation matrix BEFORE plotting or labeling
R_rounded = round(R, 2);


%% === HEATMAP (6 x 6) ===

figure('Color','w');

% Use the rounded matrix here
h = heatmap(allNames, allNames, R_rounded);

h.Title  = 'Correlation Heatmap: Spontaneous Pain vs Unpleasantness Ratings';
h.XLabel = '';
h.YLabel = '';

h.FontSize = 12;
h.GridVisible = 'on';
h.ColorbarVisible = 'on';
h.ColorLimits = [-1 1];
h.Colormap = parula;

% === CUSTOM CELL LABELS WITH BOLD FOR SIGNIFICANT CORRELATIONS ===

% Build label strings
n = size(R_rounded, 1);
labels = cell(n, n);
for r = 1:n
    for c = 1:n
        val = R_rounded(r,c);
        if isnan(val)
            labels{r,c} = '';
        else
            labels{r,c} = sprintf('%.2f', val);
        end
    end
end

% Turn off default numeric labels, overlay custom text
h.CellLabelColor = 'none';

% Access the underlying axes of the heatmap
sh = struct(h);
ax = sh.Axes;

for r = 1:n
    for c = 1:n
        if isempty(labels{r,c}), continue; end

        if P(r,c) < 0.05
            fw = 'bold';
        else
            fw = 'normal';
        end

        text(ax, c, r, labels{r,c}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontWeight', fw, ...
            'FontSize', 11, ...
            'Color','k');
    end
end


end 