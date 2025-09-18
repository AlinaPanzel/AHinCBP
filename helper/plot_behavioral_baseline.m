function plot_behavioral_baseline(d)


% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Auditory");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [0 100];

csl = [d.G1.S1.acute_mean_sound_lo;d.G2.S1.acute_mean_sound_lo;d.G3.S1.acute_mean_sound_lo];
csh = [d.G1.S1.acute_mean_sound_hi;d.G2.S1.acute_mean_sound_hi;d.G3.S1.acute_mean_sound_hi];
hsl = [d.HC.S1.acute_mean_sound_lo];
hsh = [d.HC.S1.acute_mean_sound_hi];

cpl = [d.G1.S1.acute_mean_thumb_lo;d.G2.S1.acute_mean_thumb_lo;d.G3.S1.acute_mean_thumb_lo];
cph = [d.G1.S1.acute_mean_thumb_hi;d.G2.S1.acute_mean_thumb_hi;d.G3.S1.acute_mean_thumb_hi];
hpl = [d.HC.S1.acute_mean_thumb_lo];
hph = [d.HC.S1.acute_mean_thumb_hi];

hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window


%nexttile(subplotPosition)
hc_col = [0 0 1]; 

%[0.3255 0.6863 0.6941]cyan [0 0.4470 0.7410]blue [0 0.4470 0.7410];

cbp_col = [0.8902 0.3490 0.1569];

%[0.1216 0.2863 0.0] green  [0.0 0.0627 0.4353] [0.8627 0.5647 0.2627] yellow [0.9137 0.6824 0.3020][0.6350 0.0780 0.1840];

al_goodplot(hpl, 1, 0.5, hc_col, 'left', [], std(hpl) / 1000);
al_goodplot(cpl, 1, 0.5, cbp_col, 'right', [], std(cpl) / 1000);
al_goodplot(hph, 2, 0.5, hc_col, 'left', [], std(hph) / 1000);
al_goodplot(cph, 2, 0.5, cbp_col, 'right', [], std(cph) / 1000);

al_goodplot(hsl, 4, 0.5, hc_col, 'left', [], std(hsl) / 1000);
al_goodplot(csl, 4, 0.5, cbp_col, 'right', [], std(csl) / 1000);
al_goodplot(hsh, 5, 0.5, hc_col, 'left', [], std(hsh) / 1000);
al_goodplot(csh, 5, 0.5, cbp_col, 'right', [], std(csh) / 1000);

ylim(yLimit);

% Add 'Pressure' and 'Sound' labels within the graph
text(1.5, max(ylim) * 0.92, 'Pressure', 'FontSize', 25, 'HorizontalAlignment', 'center', 'Color', 'k');
text(4.5, max(ylim) * 0.92, 'Sound', 'FontSize', 25, 'HorizontalAlignment', 'center', 'Color', 'k');

% Add a light grey vertical line at x = 3
line([3 3], ylim, 'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 2);

xticks([1 2 4 5]);
xticklabels({'Low', 'High', 'Low', 'High'});
set(gca, 'FontSize', 25); % Increased xticklabel font size
%title(titleText, 'FontSize', 30); % Increased title font size


end 