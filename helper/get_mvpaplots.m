function get_mvpaplots(lo)


% FM  

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Auditory");
%sgt.FontSize = 38;
% Define y-limits for the plots
yLimit = [-0.08 0.1];
% Plot data for different regions
plot_figures(lo.fm_mvpa, [1 2 3 4], 'FM PAIN','mvpa', yLimit, 2, false);
plot_figures(lo.fm_mvpa, [5 6 7 8], 'FM MSS','mvpa', yLimit, 2, true);




% NA

% Set up the figure
hf = figure; % Open figure and keep handle
hf = colordef(hf, 'white'); % Set color scheme
hf.Color = 'w'; % Set background color of figure window
%sgt = sgtitle("Unisensory Auditory");
%sgt.FontSize = 38;
% Define y-limits for the plots
yLimit = [-0.3 0.9];
% Plot data for different regions
plot_figures(lo.na_mvpa, [1 2 3 4], 'General','mvpa', yLimit, 3, false);
plot_figures(lo.na_mvpa, [5 6 7 8], 'Mechanical','mvpa', yLimit, 3, true);
plot_figures(lo.na_mvpa, [9 10 11 12], 'Sound','mvpa', yLimit, 3, true);





end 