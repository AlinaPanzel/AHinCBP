function get_roiplots(lo)


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
plot_figures(lo.roi, [1 2 3 4], 'Auditory Cortex','roi', yLimit, 3, true);
plot_figures(lo.roi, [5 6 7 8], 'Medial Geniculate Nucleus','roi', yLimit, 3, true);
plot_figures(lo.roi, [9 10 11 12], 'Inferior Colliculus','roi', yLimit, 3, true);


% MECHANICAL

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Mechanical");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [-1.5 2.5];

% Plot data for different regions
plot_figures(lo.roi, [1 2 3 4], 'Left Primary Somatosensory Cortex','roi', yLimit, 2, true);
plot_figures(lo.roi, [5 6 7 8], 'Right Primary Somatosensory Cortex','roi', yLimit, 2, true);


% SENSORY-INTEGRATIVE 

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Sensory-Integrative");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [-0.8 1.8];

% Plot data for different regions
plot_figures(lo.roi, [1 2 3 4], 'Ventral-anterior Insula','roi', yLimit, 3, false);
plot_figures(lo.roi, [5 6 7 8], 'Dorsal-anterior Insula','roi', yLimit, 3, true);
plot_figures(lo.roi, [9 10 11 12], 'Posterior Insula','roi', yLimit, 3, true);


% SELF-REFERENTIAL

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Self-Referential");
%sgt.FontSize = 38;

% Define y-limits for the plots
yLimit = [-1.5 1.5];

% Plot data for different regions
plot_figures(lo.roi, [1 2 3 4], 'Medial Prefrontal Cortex','roi', yLimit, 3, false);
plot_figures(lo.roi, [5 6 7 8], 'Precuneus','roi', yLimit, 3, true);
plot_figures(lo.roi, [9 10 11 12], 'Posterior Cingulate Cortex','roi', yLimit, 3, true);

end 