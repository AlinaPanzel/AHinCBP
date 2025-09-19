function f = plot_auditory_treatmenteffects_collapsed(d)
% Collapses low+high auditory unpleasantness within each session (per subject),
% then plots group means with 95% CI ribbons and SEM error bars from
% Baseline (S1) to Post (S2).
% Groups in struct d: G1=PRT, G2=Placebo, G3=UC. Each G?.S1 / G?.S2 tables
% must have columns: acute_mean_sound_lo, acute_mean_sound_hi.

alpha = 0.05; % for 95% CI

% ---- per-group/session subject-wise averages (collapse lo+hi) ----
[mPRT,ciPRT,sePRT] = grp_session_stats(d.G1, alpha);
[mPLA,ciPLA,sePLA] = grp_session_stats(d.G2, alpha);
[mUC, ciUC, seUC]  = grp_session_stats(d.G3, alpha);

% ---- build structs explicitly (avoids dot-indexing into non-structs) ----
m  = struct('PRT',     [mPRT.S1,     mPRT.S2], ...
            'placebo', [mPLA.S1,     mPLA.S2], ...
            'UC',      [mUC.S1,      mUC.S2]);
ci = struct('PRT',     [ciPRT.S1,    ciPRT.S2], ...
            'placebo', [ciPLA.S1,    ciPLA.S2], ...
            'UC',      [ciUC.S1,     ciUC.S2]);
se = struct('PRT',     [sePRT.S1,    sePRT.S2], ...
            'placebo', [sePLA.S1,    sePLA.S2], ...
            'UC',      [seUC.S1,     seUC.S2]);

% ---- plot ----
x = [1 2];                                % Baseline, Post-Tx
col.PRT     = [0.25 0.55 0.95];           % blue
col.placebo = [0.20 0.75 0.35];           % green
col.UC      = [0.95 0.35 0.35];           % red
offs      = [-0.08, 0, 0.08];             % jitter per group
alphaFill = 0.050; lwLine=2.6; lwErr=1.8; capSz=14; ms=7.5; mEdge=1.6;

% figure/axes
f = figure('Color','w','Position',[100 100 760 520], 'Tag','auditory_treatmenteffects_collapsed');
ax = axes('Parent',f); hold(ax,'on');

% ribbons = CI, errorbars = SEM
h1 = drawGroup(ax, x, m.PRT,     ci.PRT,     se.PRT,     col.PRT,     offs(1), alphaFill, lwLine, lwErr, capSz, ms, mEdge);
h2 = drawGroup(ax, x, m.placebo, ci.placebo, se.placebo, col.placebo, offs(2), alphaFill, lwLine, lwErr, capSz, ms, mEdge);
h3 = drawGroup(ax, x, m.UC,      ci.UC,      se.UC,      col.UC,      offs(3), alphaFill, lwLine, lwErr, capSz, ms, mEdge);

% y-lims/ticks
ylim(ax, [30 65]); yticks(ax, 30:5:65);

% main axes: black labels + light grid (no box)
xlim(ax, [0.75 2.25]);
set(ax,'XTick',x,'XTickLabel',{'Baseline','Post-Tx'});
set(ax,'FontName','Helvetica','FontSize',16,'LineWidth',1.2, ...
        'Box','off','XGrid','on','YGrid','on', ...
        'GridColor',[0.88 0.88 0.88],'GridAlpha',0.55, ...
        'XColor','k','YColor','k');
ylabel(ax,'Unpleasantness Ratings (VAS 0–100)','Color','k');

% overlay a light-gray frame only
axFrame = axes('Position',ax.Position,'Color','none', ...
               'XLim',ax.XLim,'YLim',ax.YLim, ...
               'XTick',[],'YTick',[], 'Box','on', ...
               'XColor',[0.88 0.88 0.88],'YColor',[0.88 0.88 0.88], ...
               'LineWidth',1.2,'HitTest','off','HandleVisibility','off');
uistack(ax,'top');

% stash numbers
f.UserData.means = m;
f.UserData.CI95  = ci;
f.UserData.SEM   = se;

% ---------------- local helpers ----------------
    function [mOut, ciOut, seOut] = grp_session_stats(G, a)
        % Compute per-subject average of lo/hi within each session, then group mean,
        % 95% CI half-width, and SEM.
        mOut  = struct('S1',NaN,'S2',NaN);
        ciOut = struct('S1',NaN,'S2',NaN);
        seOut = struct('S1',NaN,'S2',NaN);
        for ses = ["S1","S2"]
            if ~isfield(G, ses), continue; end
            S = G.(ses);
            if ~istable(S) || ~all(ismember({'acute_mean_sound_lo','acute_mean_sound_hi'}, S.Properties.VariableNames))
                continue
            end
            subjVal = mean([S.acute_mean_sound_lo(:), S.acute_mean_sound_hi(:)], 2, 'omitnan');
            subjVal = subjVal(~isnan(subjVal));
            n  = numel(subjVal);
            if n==0, continue; end
            mu = mean(subjVal);
            se = std(subjVal)/sqrt(n);
            hw = tinv(1-a/2, max(n-1,1)) * se;  % 95% CI half-width
            mOut.(ses)  = mu;
            ciOut.(ses) = hw;
            seOut.(ses) = se;
        end
    end

    function hLine = drawGroup(ax, x, m, ci, se, col, off, aFill, lwLine, lwErr, capSz, ms, mEdge)
        xj = x + off;
        % CI ribbon
        xu = [xj, fliplr(xj)];
        yu = [m + ci, fliplr(m - ci)];
        patch('XData',xu,'YData',yu,'FaceColor',col,'FaceAlpha',aFill, ...
              'EdgeColor','none','Parent',ax);
        % mean line + SEM error bars + markers
        hLine = plot(ax, xj, m, '-', 'Color', col, 'LineWidth', lwLine);
        errorbar(ax, xj, m, se, 'LineStyle','none', 'Color', col, ...
                 'LineWidth', lwErr, 'CapSize', capSz);
        plot(ax, xj, m, 'o', 'MarkerSize', ms, ...
             'MarkerFaceColor','w','MarkerEdgeColor',col,'LineWidth',mEdge);
    end
end
