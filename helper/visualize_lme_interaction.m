function visualize_lme_interaction(T, measureName)
% Visualize LMM with random time slopes for a single measure
% - Observed means ± SEM by group × time (intensity: dashed=low, solid=high)
% - Population-level model predictions (thick lines)
% - Caterpillar plot of subject-specific random slopes (timepoint)

    % -------- slice & typing --------
    if ~iscategorical(T.measure), T.measure = categorical(T.measure); end
    Tk = T(T.measure == categorical(string(measureName)), :);
    if isempty(Tk)
        warning('No rows for measure "%s".', measureName); return;
    end

    Tk.subID     = categorical(Tk.subID);
    Tk.group     = categorical(Tk.group);   % keep group coding
    Tk.timepoint = double(Tk.timepoint);    % 1=baseline, 2=post
    Tk.intensity = double(Tk.intensity);    % 1=low, 2=high
    Tk.value     = double(Tk.value);

    % -------- fit LMM (random slope for timepoint) --------
    lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)');

    % -------- observed summaries via GRPSTATS (mean/std/n) --------
    % grpstats returns variable names like mean_value, std_value, numel_value
    S = grpstats(Tk, {'group','timepoint','intensity'}, {'mean','std','numel'}, 'DataVars','value');
    S.sem = S.std_value ./ sqrt(S.numel_value);

    grpLevels = categories(Tk.group);
    C = lines(numel(grpLevels));

    % -------- model predictions (marginal over intensity) --------
    [Gg, Gt, Gi] = ndgrid(grpLevels, [1 2], [1 2]);   % group × time × intensity
    newT = table( categorical(Gg(:), grpLevels),  double(Gt(:)),  double(Gi(:)), ...
                  'VariableNames', {'group','timepoint','intensity'});
    % provide a dummy subID; Conditional=false => population-level
    newT.subID = categorical(repmat("dummy",height(newT),1));
    yhat = predict(lme, newT, 'Conditional', false);

    Pred = newT; Pred.yhat = yhat;
    % average over intensity to get marginal predictions
    Smod = grpstats(Pred, {'group','timepoint'}, {'mean'}, 'DataVars','yhat');
    Smod.Properties.VariableNames{'mean_yhat'} = 'yhat_mean';

    % -------- figure layout --------
    figure('Color','w','Name',sprintf('%s — LMM interaction', measureName));
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % ===== Left: Observed ± SEM with model fits =====
    ax1 = nexttile; hold(ax1,'on'); box(ax1,'off'); grid(ax1,'on');

    for gi = 1:numel(grpLevels)
        % intensity 1 (low) dashed; 2 (high) solid
        for ii = 1:2
            rows = S.group==grpLevels{gi} & S.intensity==ii;
            if ~any(rows), continue; end
            x = S.timepoint(rows);
            [x,ord] = sort(x);
            y = S.mean_value(rows); y = y(ord);
            e = S.sem(rows);        e = e(ord);

            ls = '--'; lw = 1.5; if ii==2, ls='-'; lw=2; end
            errorbar(ax1, x + (ii-1.5)*0.04, y, e, ...
                'LineStyle', ls, 'LineWidth', lw, 'Color', C(gi,:), ...
                'Marker','o', 'MarkerFaceColor', C(gi,:));
        end

        % overlay model marginal prediction
        r2 = Smod.group==grpLevels{gi};
        x2 = Smod.timepoint(r2);
        [x2,ord2] = sort(x2);
        y2 = Smod.yhat_mean(r2); y2 = y2(ord2);
        plot(ax1, x2, y2, '-', 'Color', C(gi,:), 'LineWidth', 3);
    end

    xticks(ax1,[1 2]); xticklabels(ax1, {'Baseline','Post'});
    xlabel(ax1,'Time'); ylabel(ax1, sprintf('%s (a.u.)', measureName));
    title(ax1, sprintf('%s — Observed ± SEM & LMM fit', measureName));
    legend(ax1, arrayfun(@(g)sprintf('Group %s', string(g)), grpLevels, 'UniformOutput', false), ...
        'Location','bestoutside');

    % ===== Right: random slopes (timepoint) =====
    ax2 = nexttile; box(ax2,'off'); grid(ax2,'on'); hold(ax2,'on');

    REraw = randomEffects(lme);   % table on newer MATLAB, dataset on older
    RE = [];
    % Try several conversions without relying on dataset2table
    if istable(REraw)
        RE = REraw;
    else
        try
            % Attempt generic conversion paths that often work
            RE = table(REraw);                 % sometimes works
        catch
            try
                RE = struct2table(struct(REraw));  % fallback
            catch
                RE = []; % give up, we’ll skip plotting
            end
        end
    end

    if ~isempty(RE) && istable(RE)
        % Identify the column that stores the term name
        nameCol = intersect({'Name','Term'}, RE.Properties.VariableNames, 'stable');
        if isempty(nameCol)
            % try any var that contains 'name' (case-insensitive)
            guess = find(contains(lower(RE.Properties.VariableNames),'name'),1);
            if ~isempty(guess), nameCol = RE.Properties.VariableNames(guess); end
        end

        % Identify the column that stores the estimate
        estCol  = intersect({'Estimate','est','Value'}, RE.Properties.VariableNames, 'stable');
        if isempty(estCol)
            guess = find(contains(lower(RE.Properties.VariableNames),'est'),1);
            if ~isempty(guess), estCol = RE.Properties.VariableNames(guess); end
        end

        % Find rows corresponding to the random slope on 'timepoint'
        isTimeSlope = false(height(RE),1);
        if ~isempty(nameCol)
            try
                isTimeSlope = strcmp(string(RE.(nameCol{1})),'timepoint');
            catch
                % tolerate types without error
            end
        end

        if any(isTimeSlope) && ~isempty(estCol)
            est = double(RE.(estCol{1})(isTimeSlope));
            [estS, idx] = sort(est, 'ascend');
            stem(ax2, estS, 'filled');
            yline(ax2, 0, '--', 'Color',[.6 .6 .6]);
            xlabel(ax2, 'Subjects (sorted)'); ylabel(ax2, 'Random slope: timepoint');
            title(ax2, 'Subject-specific slopes (BL→Post)');
            xlim(ax2, [0.5 numel(estS)+0.5]);
        else
            text(0.5,0.5,'Could not parse random timepoint slopes','Units','normalized', ...
                'HorizontalAlignment','center'); axis(ax2,'off');
        end
    else
        text(0.5,0.5,'Random-effects table unavailable','Units','normalized', ...
            'HorizontalAlignment','center'); axis(ax2,'off');
    end

    title(tl, sprintf('LMM (1 + timepoint | subID) — %s', measureName));
end
