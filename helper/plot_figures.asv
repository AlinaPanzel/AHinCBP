function plot_figures(regionData, subplotPosition, titleText, type, yLimit, nplot, hideYLabels)

if strcmp(type, 'behvavioral')

    csl = rmoutliers([regionData.G1.S1.acute_mean_sound_lo;regionData.G2.S1.acute_mean_sound_lo;regionData.G3.S1.acute_mean_sound_lo], 'mean'); 
    csh = rmoutliers([regionData.G1.S1.acute_mean_sound_hi;regionData.G2.S1.acute_mean_sound_hi;regionData.G3.S1.acute_mean_sound_hi], 'mean'); 
    hsl = rmoutliers([regionData.HC.S1.acute_mean_sound_lo], 'mean');
    hsh = rmoutliers([regionData.HC.S1.acute_mean_sound_hi], 'mean'); 
    
    cpl = rmoutliers([regionData.G1.S1.acute_mean_thumb_lo;regionData.G2.S1.acute_mean_thumb_lo;regionData.G3.S1.acute_mean_thumb_lo], 'mean'); 
    cph = rmoutliers([regionData.G1.S1.acute_mean_thumb_hi;regionData.G2.S1.acute_mean_thumb_hi;regionData.G3.S1.acute_mean_thumb_hi], 'mean'); 
    hpl = rmoutliers([regionData.HC.S1.acute_mean_thumb_lo], 'mean'); 
    hph = rmoutliers([regionData.HC.S1.acute_mean_thumb_hi], 'mean'); 
    
    sgt = sgtitle("Behavioural Unpleasantness Ratings"); 
    sgt.FontSize = 38; % Define y-limits for the plots yLimit = [0 100];


elseif strcmp(type, 'roi')

    % Primary sensory
    if strcmp(titleText, 'Auditory Cortex')
        subfield = "A1";
    elseif strcmp(titleText,'Medial Geniculate Nucleus')
        subfield = "MGN";
    elseif strcmp(titleText, 'Inferior Colliculus')
        subfield = "IC";

    elseif strcmp(titleText, 'Left Primary Somatosensory Cortex')
        subfield = "S1_L";
    elseif strcmp(titleText, 'Right Primary Somatosensory Cortex')
        subfield = "S1_R";

    % Sensory integrative
    elseif strcmp(titleText,  'Ventral-anterior Insula')
        subfield = "m_ventral_insula";
    elseif strcmp(titleText, 'Dorsal-anterior Insula')
        subfield = "m_dorsal_insula";
    elseif strcmp(titleText, 'Posterior Insula')
        subfield = "m_posterior_insula";

    % Self referential
    elseif strcmp(titleText, 'Medial Prefrontal Cortex')
        subfield = "mPFC";
    elseif strcmp(titleText, 'Precuneus')
        subfield = "precuneus";
    elseif strcmp(titleText, 'Posterior Cingulate Cortex')
        subfield = "PCC";

    end

    csl = rmoutliers([regionData.G1.S1.(subfield).all_s_l; regionData.G2.S1.(subfield).all_s_l; regionData.G3.S1.(subfield).all_s_l], 'mean'); 
    csh = rmoutliers([regionData.G1.S1.(subfield).all_s_h; regionData.G2.S1.(subfield).all_s_h; regionData.G3.S1.(subfield).all_s_h], 'mean'); 
    hsl = rmoutliers([regionData.HC.S1.(subfield).all_s_l], 'mean'); 
    hsh = rmoutliers([regionData.HC.S1.(subfield).all_s_h], 'mean'); 
    cpl = rmoutliers([regionData.G1.S1.(subfield).all_t_l; regionData.G2.S1.(subfield).all_t_l; regionData.G3.S1.(subfield).all_t_l], 'mean'); 
    cph = rmoutliers([regionData.G1.S1.(subfield).all_t_h; regionData.G2.S1.(subfield).all_t_h; regionData.G3.S1.(subfield).all_t_h], 'mean'); 
    hpl = rmoutliers([regionData.HC.S1.(subfield).all_t_l], 'mean'); 
    hph = rmoutliers([regionData.HC.S1.(subfield).all_t_h], 'mean'); 


elseif strcmp(type, 'mvpa')

    if strcmp(titleText, 'FM PAIN')
        subfield = "FM_PAIN"; 
    elseif strcmp(titleText, 'FM MSS')
        subfield = "FM_MSS";
    elseif strcmp(titleText, 'General')
        subfield = "general";
    elseif strcmp(titleText, 'Mechanical')
        subfield = "mechanical";
    elseif strcmp(titleText, 'Sound')
        subfield = "sound";
    end 

    csl = rmoutliers([regionData.G1.S1.all_s_l.(subfield); regionData.G2.S1.all_s_l.(subfield); regionData.G3.S1.all_s_l.(subfield)], 'mean'); 
    csh = rmoutliers([regionData.G1.S1.all_s_l.(subfield); regionData.G2.S1.all_s_l.(subfield); regionData.G3.S1.all_s_l.(subfield)], 'mean'); 
    hsl = rmoutliers([regionData.HC.S1.all_s_l.(subfield)], 'mean'); 
    hsh = rmoutliers([regionData.HC.S1.all_s_l.(subfield)], 'mean'); 

    cpl = rmoutliers([regionData.G1.S1.all_t_l.(subfield); regionData.G2.S1.all_t_l.(subfield); regionData.G3.S1.all_t_l.(subfield)], 'mean'); 
    cph = rmoutliers([regionData.G1.S1.all_t_l.(subfield); regionData.G2.S1.all_t_l.(subfield); regionData.G3.S1.all_t_l.(subfield)], 'mean'); 
    hpl = rmoutliers([regionData.HC.S1.all_t_l.(subfield)], 'mean'); 
    hph = rmoutliers([regionData.HC.S1.all_t_l.(subfield)], 'mean'); 

end


% Determine the number of subplots
if nplot == 2
    nsubplot = 8;
elseif nplot == 3
    nsubplot = 12;
end
% Set up tiled layout
%t = tiledlayout()


% Plot data
if strcmp(type, 'behavioral')
    disp('no subplot')
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
% 
 else 
 subplot(1, nsubplot, subplotPosition)
 end

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
set(gca, 'FontSize', 18); % Increased xticklabel font size
title(titleText, 'FontSize', 27); % Increased title font size

% Adjust y-axis labels if needed
%if hideYLabels
 %   yticklabels([]); % Hide y-axis labels
%else 
%    ax = gca;
%labels = string(ax.YAxis.TickLabels); % extract
%labels(2:2:end) = nan; % remove every other one
%ax.YAxis.TickLabels = labels; % set
end
