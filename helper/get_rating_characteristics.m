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



% Extract variables from the summaryTable
sp_mean = summaryTable.spon_mean;
sp_var  = summaryTable.spon_var;

sound_low    = summaryTable.sound_low_mean;
sound_high   = summaryTable.sound_high_mean;
press_low    = summaryTable.press_low_mean;
press_high   = summaryTable.press_high_mean;

% Compute correlations (Pearson)
[r_spmean_sound_low,    p_spmean_sound_low]    = corr(sp_mean, sound_low,  'rows','pairwise');
[r_spmean_sound_high,   p_spmean_sound_high]   = corr(sp_mean, sound_high, 'rows','pairwise');
[r_spmean_press_low,    p_spmean_press_low]    = corr(sp_mean, press_low,  'rows','pairwise');
[r_spmean_press_high,   p_spmean_press_high]   = corr(sp_mean, press_high, 'rows','pairwise');

[r_spvar_sound_low,     p_spvar_sound_low]     = corr(sp_var, sound_low,  'rows','pairwise');
[r_spvar_sound_high,    p_spvar_sound_high]    = corr(sp_var, sound_high, 'rows','pairwise');
[r_spvar_press_low,     p_spvar_press_low]     = corr(sp_var, press_low,  'rows','pairwise');
[r_spvar_press_high,    p_spvar_press_high]    = corr(sp_var, press_high, 'rows','pairwise');

% Print results to command window
fprintf('\n=== CORRELATIONS (Reviewer I) ===\n');
fprintf('Spont pain MEAN  vs Sound Low :   r = %.2f, p = %.3f\n', r_spmean_sound_low,  p_spmean_sound_low);
fprintf('Spont pain MEAN  vs Sound High:   r = %.2f, p = %.3f\n', r_spmean_sound_high, p_spmean_sound_high);
fprintf('Spont pain MEAN  vs Press Low :   r = %.2f, p = %.3f\n', r_spmean_press_low,  p_spmean_press_low);
fprintf('Spont pain MEAN  vs Press High:   r = %.2f, p = %.3f\n', r_spmean_press_high, p_spmean_press_high);

fprintf('Spont pain VAR   vs Sound Low :   r = %.2f, p = %.3f\n', r_spvar_sound_low,  p_spvar_sound_low);
fprintf('Spont pain VAR   vs Sound High:   r = %.2f, p = %.3f\n', r_spvar_sound_high, p_spvar_sound_high);
fprintf('Spont pain VAR   vs Press Low :   r = %.2f, p = %.3f\n', r_spvar_press_low,  p_spvar_press_low);
fprintf('Spont pain VAR   vs Press High:   r = %.2f, p = %.3f\n', r_spvar_press_high, p_spvar_press_high);




%% === HEATMAP (6 x 6) ===

 % Round full correlation matrix BEFORE plotting or labeling
R_rounded = round(R, 2);

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

% Significance stars added on top of already rounded data
labels = strings(size(R_rounded));
for r = 1:size(R_rounded,1)
    for c = 1:size(R_rounded,2)
        if P(r,c) < 0.05
            labels(r,c) = sprintf('%.2f*', R_rounded(r,c));
        else
            labels(r,c) = sprintf('%.2f', R_rounded(r,c));
        end
    end
end

h.CellLabelFormat = '%s';
h.CellLabelData   = labels;

end 