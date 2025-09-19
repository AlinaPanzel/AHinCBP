function f = plot_mpfc_treatmenteffects_sound(T)
% Plot mPFC treatment effects (Sound only), collapsing Low/High per session.
% Ribbons = 95% CI; Error bars = SEM.
% T: columns subID, timepoint(1/2), group(1=PRT,2=Placebo,3=UC), intensity, measure, value

% --- filter to mPFC (sound only table in your setup) ---
T.measure = string(T.measure);
Tk = T(T.measure=="mPFC",:);

% --- collapse Low/High within subject x session ---
G = groupsummary(Tk, {'subID','group','timepoint'}, 'mean', 'value');  % mean across intensities

% --- compute means, 95% CI (half-width), and SEM per group/session ---
alpha  = 0.05;                          % 95% CI
grpVals = [1 2 3];                      % 1=PRT, 2=Placebo, 3=UC
labels  = {'PRT','Placebo','UC'};
means   = nan(3,2);                     % rows=group, cols=S1,S2
cihw    = nan(3,2);                     % 95% CI half-width (for ribbon)
sems    = nan(3,2);                     % SEM (for error bars)

for gi = 1:numel(grpVals)
    g = grpVals(gi);
    for s = 1:2
        v = G.mean_value(G.group==g & G.timepoint==s);
        v = v(~isnan(v));
        n = numel(v);
        if n>0
            mu = mean(v);
            se = std(v)/sqrt(n);
            hw = tinv(1-alpha/2, max(n-1,1)) * se;   % 95% CI half-width
            means(gi,s) = mu;
            cihw(gi,s)  = hw;
            sems(gi,s)  = se;
        end
    end
end

% --- plot ---
x      = [1 2];                 % Baseline, Post
offs   = [-0.08 0 0.08];        % jitter per group
cols   = [0.25 0.55 0.95; 0.20 0.75 0.35; 0.95 0.35 0.35];   % PRT/Placebo/UC

f  = figure('Color','w','Position',[120 120 760 520], 'Tag','mPFC_treatmenteffects_sound');
ax = axes('Parent',f); hold(ax,'on');

lwLine=2.6; lwErr=1.8; capSz=14; ms=7.5; mEdge=1.6; alphaFill=0.05;
h = gobjects(1,3);

for gi = 1:3
    xj = x + offs(gi);
    m  = means(gi,:);
    hw = cihw(gi,:);   % 95% CI for ribbons
    se = sems(gi,:);   % SEM for error bars
    c  = cols(gi,:);

    % 95% CI ribbon
    xu = [xj, fliplr(xj)];
    yu = [m + hw, fliplr(m - hw)];
    patch('XData',xu,'YData',yu,'FaceColor',c,'FaceAlpha',alphaFill, ...
          'EdgeColor','none','Parent',ax);

    % mean line + SEM error bars + markers
    h(gi) = plot(ax, xj, m, '-', 'Color', c, 'LineWidth', lwLine);
    errorbar(ax, xj, m, se, 'LineStyle','none', 'Color', c, ...
             'LineWidth', lwErr, 'CapSize', capSz);
    plot(ax, xj, m, 'o', 'MarkerSize', ms, 'MarkerFaceColor','w', ...
         'MarkerEdgeColor', c, 'LineWidth', mEdge);
end

% ----- axes cosmetics: light grid + black labels, y starts at -0.25 -----
xlim(ax,[0.75 2.25]); set(ax,'XTick',[1 2],'XTickLabel',{'Baseline','Post-Tx'});
yl = [min(means-cihw,[],'all','omitnan') max(means+cihw,[],'all','omitnan')];
if ~all(isfinite(yl)), yl = [-0.25 0.06]; end
yl(1) = -0.25; pad = 0.05*range(yl);
ylim(ax,[yl(1) yl(2)+pad]);

set(ax,'FontName','Helvetica','FontSize',16,'LineWidth',1.2, ...
        'Box','off','XGrid','on','YGrid','on', ...
        'GridColor',[0.88 0.88 0.88],'GridAlpha',0.55, ...
        'XColor','k','YColor','k');
ylabel(ax,'mPFC activity','Color','k');

% light gray frame box on top
axFrame = axes('Position',ax.Position,'Color','none','XLim',ax.XLim,'YLim',ax.YLim, ...
               'XTick',[],'YTick',[],'Box','on','XColor',[0.85 0.85 0.85], ...
               'YColor',[0.85 0.85 0.85],'LineWidth',1.2,'HitTest','off','HandleVisibility','off');
uistack(ax,'top');

% legend if you want it here:
% legend(ax,h,labels,'Orientation','horizontal','Location','northoutside','Box','off');

% stash numbers
f.UserData.means = struct('PRT',means(1,:), 'Placebo',means(2,:), 'UC',means(3,:));
f.UserData.CI95  = struct('PRT',cihw(1,:),  'Placebo',cihw(2,:),  'UC',cihw(3,:));
f.UserData.SEM   = struct('PRT',sems(1,:),  'Placebo',sems(2,:),  'UC',sems(3,:));
end
