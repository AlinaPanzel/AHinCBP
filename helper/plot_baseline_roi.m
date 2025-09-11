function plot_baseline_roi

% AUDITORY 

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Auditory");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [-2 7];

% Plot data for different regions
plot_allgoodplot(roi.A1, [1 2 3 4], 'Auditory Cortex','roi', yLimit, 3, true);
plot_allgoodplot(roi.MGN, [5 6 7 8], 'Medial Geniculate Nucleus','roi', yLimit, 3, true);
plot_allgoodplot(roi.IC, [9 10 11 12], 'Inferior Colliculus','roi', yLimit, 3, true);


end 
