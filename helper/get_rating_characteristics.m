function baseline_ratings = get_rating_characteristics(d)

ids   = d.metadata.id;
spon  = d.metadata.spon_pain_ratings;      % N x 8
unpl  = d.metadata.acute_pain_ratings;     % N x 10
thumb = d.metadata.acute_stim_is_thumb;    % N x 10 (0=sound, 1=thumb)
inten = d.metadata.acute_stim_intensity;   % N x 10 (1=low, 2=high)

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

figure;

% Use the rounded matrix here
h = heatmap(allNames, allNames, R_rounded);

h.Title  = 'Correlation Heatmap: Spontaneous & Stimulus Ratings';
h.XLabel = 'Variables';
h.YLabel = 'Variables';

h.FontSize = 12;
h.GridVisible = 'off';
h.ColorbarVisible = 'on';
h.ColorLimits = [-1 1];

% Alternative color palette
try
    h.Colormap = cividis;    % colorblind-friendly, perceptually uniform
catch
    h.Colormap = parula;
end

set(gcf,'Color','w');
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

% === CUSTOM CELL LABELS WITH BOLD FOR SIGNIFICANT CORRELATIONS ===

% Turn off built-in labels (we'll draw them manually)
h.CellLabelColor = 'none';

% Access the underlying axes of the heatmap
sh = struct(h);
ax = sh.Axes;

n = size(R_rounded, 1);

for r = 1:n
    for c = 1:n

        val = R_rounded(r,c);
        if isnan(val)
            continue; % skip NaNs
        end

        if P(r,c) < 0.05
            txt = sprintf('%.2f*', val);
            fw  = 'bold';    % significant = bold
        else
            txt = sprintf('%.2f', val);
            fw  = 'normal';  % non-significant = normal
        end

        % Draw text in the heatmap axes
        text(ax, c, r, txt, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontWeight', fw, ...
            'FontSize', 11, ...
            'Color','k');
    end
end


h.CellLabelFormat = '%s';
h.CellLabelData   = labels;


end 