function plot_auditory_treatmenteffects_byIntensity(d)
% Auditory unpleasantness by intensity and group over time.
% G1=PRT (green), G2=Placebo (red), G3=UC (blue)
% High: solid bold; Low: dashed thin. Baseline=S1 at x=1, Post=S2 at x=2.

alpha = 0.05;

% ---- summaries (means ± 95% CI) per intensity and session ----
PRT = summarize_grp(d.G1, alpha);   % struct with .lo = [S1 S2], .hi = [S1 S2], same for .ci
PLA = summarize_grp(d.G2, alpha);
UC  = summarize_grp(d.G3, alpha);

% ---- colors (locked to groups) ----
col.PRT = [0.25 0.55 0.95];   % blue
col.PLA = [0.20 0.75 0.35];   % green
col.UC  = [0.95 0.35 0.35];   % red

% ---- plot setup ----
x = [1 2];                         % Baseline, Post
offs = [-0.08, 0, 0.08];           % jitter per group
alphaFill = 0.00; lwBold=2.8; lwThin=1.8; lwErr=1.6; capSz=14; ms=7.5; mEdge=1.6;

figure('Color','w','Position',[120 100 820 540]); hold on
draw = @(m,ci,c,off,isHigh) drawLine(x,m,ci,c,off,alphaFill,lwBold,lwThin,lwErr,capSz,ms,mEdge,isHigh);

% PRT (G1)
h_prt_hi = draw(PRT.hi, PRT.ci.hi, col.PRT, offs(1), true);
h_prt_lo = draw(PRT.lo, PRT.ci.lo, col.PRT, offs(1), false);

% Placebo (G2)
h_pla_hi = draw(PLA.hi, PLA.ci.hi, col.PLA, offs(2), true);
h_pla_lo = draw(PLA.lo, PLA.ci.lo, col.PLA, offs(2), false);

% Usual Care (G3)
h_uc_hi  = draw(UC.hi,  UC.ci.hi,  col.UC,  offs(3), true);
h_uc_lo  = draw(UC.lo,  UC.ci.lo,  col.UC,  offs(3), false);

% ---- axes / labels ----
% dynamic y-lims from data
allMeans = [PRT.lo, PRT.hi, PLA.lo, PLA.hi, UC.lo, UC.hi];
allCIs   = [PRT.ci.lo, PRT.ci.hi, PLA.ci.lo, PLA.ci.hi, UC.ci.lo, UC.ci.hi];
yl = [min(allMeans - allCIs), max(allMeans + allCIs)];
if any(isfinite(yl))
    pad = 0.08*range(yl); if pad==0, pad=0.5; end
    ylim([yl(1)-pad, yl(2)+pad]);
end
xlim([0.75 2.25]);

set(gca,'XTick',x,'XTickLabel',{'Baseline','Post-Tx'}, ...
        'FontName','Helvetica','FontSize',12, 'LineWidth',1.2, 'Box','off')
ylabel('Auditory unpleasantness (VAS 0–100)','FontWeight','bold')
%grid on

% Title and legend
sgtitle('Effects of Treatment on Auditory Unpleasantness Ratings (by Intensity)', ...
        'FontSize', 16, 'FontWeight', 'bold');

lg = legend([h_prt_hi h_prt_lo h_pla_hi h_pla_lo h_uc_hi h_uc_lo], ...
            {'PRT High','PRT Low','Placebo High','Placebo Low','UC High','UC Low'}, ...
            'Orientation','horizontal', 'Location','northoutside', 'Box','off');

% ---- sanity printout so you can verify values ----
fprintf('\n== Means being plotted (Baseline=S1, Post=S2) ==\n');
lab = {'PRT','PLA','UC'};
dat = {PRT, PLA, UC};
for i = 1:3
    S = dat{i};
    fprintf('%s Low:  S1=%6.2f  S2=%6.2f   | High: S1=%6.2f  S2=%6.2f\n', ...
        lab{i}, S.lo(1), S.lo(2), S.hi(1), S.hi(2));
end
fprintf('===============================================\n');

% ----------------- helpers -----------------
function S = summarize_grp(G, a)
    % Returns S.lo, S.hi (1x2 vectors for [S1 S2])
    % and S.ci.lo, S.ci.hi (1x2 half-widths)
    S.lo = [NaN NaN]; S.hi = [NaN NaN];
    S.ci.lo = [NaN NaN]; S.ci.hi = [NaN NaN];

    % enforce S1->index1, S2->index2 mapping
    if isfield(G,'S1') && istable(G.S1)
        checkTime(G.S1, 'S1');
        S.lo(1)    = mean(G.S1.acute_mean_sound_lo,'omitnan');
        S.hi(1)    = mean(G.S1.acute_mean_sound_hi,'omitnan');
        S.ci.lo(1) = ciHalfWidth(G.S1.acute_mean_sound_lo, a);
        S.ci.hi(1) = ciHalfWidth(G.S1.acute_mean_sound_hi, a);
    end
    if isfield(G,'S2') && istable(G.S2)
        checkTime(G.S2, 'S2');
        S.lo(2)    = mean(G.S2.acute_mean_sound_lo,'omitnan');
        S.hi(2)    = mean(G.S2.acute_mean_sound_hi,'omitnan');
        S.ci.lo(2) = ciHalfWidth(G.S2.acute_mean_sound_lo, a);
        S.ci.hi(2) = ciHalfWidth(G.S2.acute_mean_sound_hi, a);
    end
end

function hw = ciHalfWidth(x, a)
    x = x(:); x = x(isfinite(x));
    n = numel(x);
    if n==0, hw = NaN; return; end
    se = std(x)/sqrt(n);
    hw = tinv(1-a/2, max(n-1,1)) * se;
end

function checkTime(T, sesName)
    % Optional sanity: warn if time column mismatches expected code
    if ismember('time', T.Properties.VariableNames)
        tt = unique(T.time);
        if any(~ismember(tt, [1 2]))
            warning('%s has time values not in {1,2}: %s', sesName, mat2str(tt));
        end
        if sesName=="S1" && any(tt~=1)
            warning('S1 contains time~=1 (found %s). Using as Baseline anyway.', mat2str(tt));
        end
        if sesName=="S2" && any(tt~=2)
            warning('S2 contains time~=2 (found %s). Using as Post anyway.', mat2str(tt));
        end
    end
end

function h = drawLine(x, m, ci, col, off, aFill, lwBold, lwThin, lwErr, capSz, ms, mEdge, isHigh)
    % style by intensity
    if isHigh
        ls = '-'; lw = lwBold;
    else
        ls = '--'; lw = lwThin;
    end
    xj = x + off;
    % shaded CI
    if all(isfinite(ci))
        patch([xj fliplr(xj)], [m+ci fliplr(m-ci)], col, 'FaceAlpha', aFill, 'EdgeColor','none');
    end
    % main line
    h = plot(xj, m, ls, 'Color', col, 'LineWidth', lw);
    % error bars
    if all(isfinite(ci))
        errorbar(xj, m, ci, 'LineStyle','none', 'Color', col, 'LineWidth', lwErr, 'CapSize', capSz);
    end
    % markers
    plot(xj, m, 'o', 'MarkerSize', ms, 'MarkerFaceColor','w', 'MarkerEdgeColor', col, 'LineWidth', mEdge);
end

end
