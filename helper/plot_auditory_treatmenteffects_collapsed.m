function f = plot_auditory_treatmenteffects_collapsed(d)

% Collapses low+high auditory unpleasantness within each session (per subject),
% then plots group means ±95% CI from Baseline (S1) to Post (S2).
% Groups: G1=PRT, G2=Placebo, G3=UC.

alpha = 0.05;

% ---- per-group/session subject-wise averages (collapse lo+hi) ----
[mPRT,ciPRT] = grp_session_meanCI(d.G1);
[mPLA,ciPLA] = grp_session_meanCI(d.G2);
[mUC, ciUC ] = grp_session_meanCI(d.G3);

% map plotting keys
m.PRT     = [mPRT.S1, mPRT.S2];     ci.PRT     = [ciPRT.S1, ciPRT.S2];
m.placebo = [mPLA.S1, mPLA.S2];     ci.placebo = [ciPLA.S1, ciPLA.S2];
m.UC      = [mUC.S1,  mUC.S2];      ci.UC      = [ciUC.S1,  ciUC.S2];

% ---- plot ----
x = [1 2]; % Baseline, Post-Tx

% colors (fixed to groups)
col.PRT     = [0.25 0.55 0.95];   % blue
col.placebo = [0.20 0.75 0.35];   % green
col.UC      = [0.95 0.35 0.35];   % red

offs      = [-0.08, 0, 0.08];  % jitter per group
alphaFill = 0.05; lwLine=2.6; lwErr=1.8; capSz=14; ms=7.5; mEdge=1.6;

% Make a fresh figure and explicit axes; return this handle
f = figure('Color','w','Position',[100 100 760 520],'Tag','auditory_treatmenteffects_collapsed');
ax = axes('Parent',f); hold(ax,'on');

plotGroup = @(meanVals,ciVals,col,off) ...
    drawGroup(ax, x, meanVals, ciVals, col, off, alphaFill, lwLine, lwErr, capSz, ms, mEdge);

h1 = plotGroup(m.PRT,     ci.PRT,     col.PRT,     offs(1));
h2 = plotGroup(m.placebo, ci.placebo, col.placebo, offs(2));
h3 = plotGroup(m.UC,      ci.UC,      col.UC,      offs(3));

% y-lims
allVals = [m.PRT, m.placebo, m.UC]; allCI = [ci.PRT, ci.placebo, ci.UC]; %#ok<NASGU>
ylim(ax, [30 65]); yticks(ax, 30:5:65);

xlim(ax, [0.75 2.25]);
set(ax,'XTick',x,'XTickLabel',{'Baseline','Post-Tx'}, ...
       'FontName','Helvetica','FontSize',12, 'LineWidth',1.2, 'Box','off')
ylabel(ax,'Auditory unpleasantness (VAS 0–100)','FontWeight','bold')

% Title above everything
sgtitle(f,'Effects of Treatment on Auditory Unpleasantness Ratings', ...
        'FontSize',16,'FontWeight','bold');

% Legend under title
lg = legend(ax, [h1 h2 h3], {'PRT','Placebo','UC'}, ...
            'Orientation','horizontal', ...
            'Location','northoutside', ...
            'Box','off'); %#ok<NASGU>

% Store summary in figure (handy when saving)
f.UserData.means = m;
f.UserData.ci    = ci;

% ---------------- local helpers ----------------
    function [mOut, ciOut] = grp_session_meanCI(G)
        % Compute per-subject average of lo/hi within each session, then mean±CI
        mOut = struct('S1',NaN,'S2',NaN);
        ciOut= struct('S1',NaN,'S2',NaN);
        for ses = ["S1","S2"]
            if ~isfield(G, ses), continue; end
            S = G.(ses);
            if ~istable(S) || ~all(ismember({'acute_mean_sound_lo','acute_mean_sound_hi'}, S.Properties.VariableNames))
                continue
            end
            subjVal = mean([S.acute_mean_sound_lo(:), S.acute_mean_sound_hi(:)], 2, 'omitnan');
            subjVal = subjVal(~isnan(subjVal));
            [mOut.(ses), ciOut.(ses)] = meanCI(subjVal, alpha);
        end
    end

    function [mu, ci] = meanCI(x, a)
        x = x(:); x = x(~isnan(x));
        n = numel(x);
        if n==0, mu=NaN; ci=NaN; return; end
        mu = mean(x);
        se = std(x)/sqrt(n);
        ci = tinv(1-a/2, max(n-1,1)) * se;  % half-width
    end

    function hLine = drawGroup(ax, x, m, ci, col, off, aFill, lwLine, lwErr, capSz, ms, mEdge)
        xj = x + off;
        xu = [xj, fliplr(xj)];
        yu = [m + ci, fliplr(m - ci)];
        patch('XData',xu,'YData',yu,'FaceColor',col,'FaceAlpha',aFill, ...
              'EdgeColor','none','Parent',ax);
        hLine = plot(ax, xj, m, '-', 'Color', col, 'LineWidth', lwLine);
        errorbar(ax, xj, m, ci, 'LineStyle','none', 'Color', col, ...
                 'LineWidth', lwErr, 'CapSize', capSz);
        plot(ax, xj, m, 'o', 'MarkerSize', ms, ...
             'MarkerFaceColor','w','MarkerEdgeColor',col,'LineWidth',mEdge);
    end
end

