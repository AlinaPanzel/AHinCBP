function hf = plot_behavioral_baseline(lo, d, type, varargin)
% PLOT_BEHAVIORAL_BASELINE(lo, d, type, ...)
% type: 'behavioural' | 'roi' | 'fm' | 'na' | 'mvpa' (mvpa == na)
%
% Examples:
%   plot_behavioral_baseline(lo, d, 'roi');        % all ROIs (tiled)
%   plot_behavioral_baseline(lo, d, 'fm');         % all FM metrics
%   plot_behavioral_baseline(lo, d, 'na');         % all NA metrics
%   plot_behavioral_baseline(lo, d, 'behavioural') % behavioural only (uses d)

% ------------------------- params -----------------------------------------
p = inputParser;
p.addOptional('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.parse(varargin{:});
yLimit = p.Results.yLim;

type = lower(string(type));

% ------------------------- figure -----------------------------------------
hf = figure('Color','w');
try, colordef(hf,'white'); end %#ok<COLND>

% Colors
hc_col  = [0 0 1];                % controls/HC
cbp_col = [0.8902 0.3490 0.1569]; % patients

% Convenience lambdas
firstField = @(S, candidates) candidates{find(cellfun(@(c) isfield(S,c), candidates), 1, 'first')};
hasField   = @(S, f) isfield(S, f);

% Groups that may exist
grpNames = {'G1','G2','G3','HC'};

switch type

%% ===================== BEHAVIOURAL: single figure =======================
case 'behavioural'
    assert(~isempty(d) && isstruct(d), 'For type="behavioural" provide non-empty struct d.');
    csl = cat_exist(d, grpNames(1:3), @(g) g.S1.acute_mean_sound_lo);
    csh = cat_exist(d, grpNames(1:3), @(g) g.S1.acute_mean_sound_hi);
    hsl = pick_exist(d, 'HC', @(g) g.S1.acute_mean_sound_lo);
    hsh = pick_exist(d, 'HC', @(g) g.S1.acute_mean_sound_hi);

    cpl = cat_exist(d, grpNames(1:3), @(g) g.S1.acute_mean_thumb_lo);
    cph = cat_exist(d, grpNames(1:3), @(g) g.S1.acute_mean_thumb_hi);
    hpl = pick_exist(d, 'HC', @(g) g.S1.acute_mean_thumb_lo);
    hph = pick_exist(d, 'HC', @(g) g.S1.acute_mean_thumb_hi);

    sgtitle('Behavioural Unpleasantness Ratings','FontSize',22);
    do_panel_plot(hc_col, cbp_col, hpl, cpl, hph, cph, hsl, csl, hsh, csh, yLimit);

%% ===================== ROI: tiled across all =======================
case 'roi'
    % locate the ROI container field name (ROI/roi/ROIs/rois) inside lo.roi.<group>.S1
    assert(isstruct(lo) && isfield(lo,'roi'), 'Expected lo.roi structure.');
    oneGroup = first_present(lo.roi, grpNames);
    assert(~isempty(oneGroup), 'No groups (G1/G2/G3/HC) found in lo.roi.');
    S1 = lo.roi.(oneGroup).S1;

    roiField = firstField(S1, {'ROI','roi','ROIs','rois'});
    assert(~isempty(roiField), 'Could not find ROI container (ROI/roi/ROIs/rois) under lo.roi.%s.S1.', oneGroup);

    nROI = numel(S1.(roiField));
    tl = tiledlayout(ceil(nROI/3), 3, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(tl, 'ROI results (all)','FontSize',22);

    for r = 1:nROI
        nexttile;

        % assemble series from groups that exist
        [csl,csh,hsl,hsh,cpl,cph,hpl,hph] = roi_series(lo, roiField, r, grpNames);

        do_panel_plot(hc_col, cbp_col, hpl, cpl, hph, cph, hsl, csl, hsh, csh, yLimit);
        title(sprintf('ROI %d', r), 'FontSize',12);
    end

%% ===================== FM: tiled across all ===================
case 'fm'
    assert(isstruct(lo) && isfield(lo,'fm_mvpa'), 'Expected lo.fm_mvpa structure.');
    % infer the available FM metric names from any present group
    root = lo.fm_mvpa;
    oneGroup = first_present(root, grpNames);
    assert(~isempty(oneGroup), 'No groups (G1/G2/G3/HC) found in lo.fm_mvpa.');

    % Find available metric fields under all_s_l.* for S1
    fmCandidates = metrics_from_struct(root.(oneGroup).S1, 'all_s_l');
    assert(~isempty(fmCandidates), 'No FM metrics found under lo.fm_mvpa.%s.S1.*', oneGroup);

    nM = numel(fmCandidates);
    tl = tiledlayout(ceil(nM/3), 3, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(tl, 'FM MVPA results (all)','FontSize',22);

    for m = 1:nM
        key = fmCandidates{m};
        nexttile;

        [csl,csh,hsl,hsh,cpl,cph,hpl,hph] = mvpa_series(root, key, grpNames);

        do_panel_plot(hc_col, cbp_col, hpl, cpl, hph, cph, hsl, csl, hsh, csh, yLimit);
        title(key, 'Interpreter','none','FontSize',12);
    end

%% ===================== NA/MVPA: tiled across all =============
case {'na','mvpa'}
    % Support NA in either lo.NA_mvpa or lo.na_mvpa (case tolerant)
    naRootField = firstField(lo, {'NA_mvpa','na_mvpa'});
    assert(~isempty(naRootField), 'Expected lo.NA_mvpa (or lo.na_mvpa) structure.');
    root = lo.(naRootField);

    oneGroup = first_present(root, grpNames);
    assert(~isempty(oneGroup), 'No groups (G1/G2/G3/HC) found in lo.%s.', naRootField);

    naCandidates = metrics_from_struct(root.(oneGroup).S1, 'all_s_l');
    assert(~isempty(naCandidates), 'No NA/MVPA metrics found under lo.%s.%s.S1.*', naRootField, oneGroup);

    nM = numel(naCandidates);
    tl = tiledlayout(ceil(nM/3), 3, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(tl, 'NA/MVPA results (all)','FontSize',22);

    for m = 1:nM
        key = naCandidates{m};
        nexttile;

        [csl,csh,hsl,hsh,cpl,cph,hpl,hph] = mvpa_series(root, key, grpNames);

        do_panel_plot(hc_col, cbp_col, hpl, cpl, hph, cph, hsl, csl, hsh, csh, yLimit);
        title(key, 'Interpreter','none','FontSize',12);
    end

otherwise
    error('type must be one of: ''behavioural'', ''roi'', ''fm'', ''na'' (or ''mvpa'').');
end

% ====================== helpers ==========================================
    function gname = first_present(rootS, names)
        gname = '';
        for ii = 1:numel(names)
            if hasField(rootS, names{ii})
                gname = names{ii};
                return;
            end
        end
    end

    function lst = metrics_from_struct(S1, prefix)
        % returns fieldnames available under S1.(prefix).*
        if ~hasField(S1, 'all_s_l') && ~hasField(S1, prefix)
            lst = {};
            return;
        end
        if hasField(S1, prefix)
            f = fieldnames(S1.(prefix));
            lst = f(:)';
        else
            f = fieldnames(S1.all_s_l);
            lst = f(:)';
        end
    end

    function arr = cat_exist(rootS, grps, getter)
        % concatenate over groups that exist (patients)
        tmp = [];
        for jj = 1:numel(grps)
            g = grps{jj};
            if hasField(rootS, g)
                try
                    tmp = [tmp; getter(rootS.(g))]; %#ok<AGROW>
                catch
                end
            end
        end
        arr = tmp;
    end

    function arr = pick_exist(rootS, gname, getter)
        if hasField(rootS, gname)
            try
                arr = getter(rootS.(gname));
            catch
                arr = [];
            end
        else
            arr = [];
        end
    end

    function [csl,csh,hsl,hsh,cpl,cph,hpl,hph] = roi_series(loS, roiFieldName, idx, grps)
        % Build series for one ROI index across available groups
        csl = []; csh = []; cpl = []; cph = [];
        for jj = 1:3 % G1..G3
            g = sprintf('G%d', jj);
            if isfield(loS.roi, g)
                R = loS.roi.(g).S1.(roiFieldName){idx};
                csl = [csl; safe_get(R,'all_s_l')]; %#ok<AGROW>
                csh = [csh; safe_get(R,'all_s_h')]; %#ok<AGROW>
                cpl = [cpl; safe_get(R,'all_p_l')]; %#ok<AGROW>
                cph = [cph; safe_get(R,'all_p_h')]; %#ok<AGROW>
            end
        end
        if isfield(loS.roi,'HC')
            R = loS.roi.HC.S1.(roiFieldName){idx};
            hsl = safe_get(R,'all_s_l');
            hsh = safe_get(R,'all_s_h');
            hpl = safe_get(R,'all_p_l');
            hph = safe_get(R,'all_p_h');
        else
            hsl=[]; hsh=[]; hpl=[]; hph=[];
        end
    end

    function [csl,csh,hsl,hsh,cpl,cph,hpl,hph] = mvpa_series(root, key, grps)
        % Build series for one metric name across available groups
        csl = cat_exist(root, grps(1:3), @(g) g.S1.all_s_l.(key));
        csh = cat_exist(root, grps(1:3), @(g) g.S1.all_s_h.(key));
        hsl = pick_exist(root, 'HC', @(g) g.S1.all_s_l.(key));
        hsh = pick_exist(root, 'HC', @(g) g.S1.all_s_h.(key));

        cpl = cat_exist(root, grps(1:3), @(g) g.S1.all_p_l.(key));
        cph = cat_exist(root, grps(1:3), @(g) g.S1.all_p_h.(key));
        hpl = pick_exist(root, 'HC', @(g) g.S1.all_p_l.(key));
        hph = pick_exist(root, 'HC', @(g) g.S1.all_p_h.(key));
    end

    function v = safe_get(S, fieldname)
        if isfield(S, fieldname)
            v = S.(fieldname);
        else
            v = [];
        end
    end

    function do_panel_plot(hc_col, cbp_col, hpl, cpl, hph, cph, hsl, csl, hsh, csh, yLimitLocal)
        % main 5-position panel with your al_goodplot calls

        % Pressure
        if ~isempty(hpl), al_goodplot(hpl, 1, 0.5, hc_col,  'left',  [], std(hpl)/1000); end
        if ~isempty(cpl), al_goodplot(cpl, 1, 0.5, cbp_col, 'right', [], std(cpl)/1000); end
        if ~isempty(hph), al_goodplot(hph, 2, 0.5, hc_col,  'left',  [], std(hph)/1000); end
        if ~isempty(cph), al_goodplot(cph, 2, 0.5, cbp_col, 'right', [], std(cph)/1000); end

        % Sound
        if ~isempty(hsl), al_goodplot(hsl, 4, 0.5, hc_col,  'left',  [], std(hsl)/1000); end
        if ~isempty(csl), al_goodplot(csl, 4, 0.5, cbp_col, 'right', [], std(csl)/1000); end
        if ~isempty(hsh), al_goodplot(hsh, 5, 0.5, hc_col,  'left',  [], std(hsh)/1000); end
        if ~isempty(csh), al_goodplot(csh, 5, 0.5, cbp_col, 'right', [], std(csh)/1000); end

        % y-limits (auto if none provided)
        if isempty(yLimitLocal)
            % Auto span based on data, with 5% padding
            allv = [hpl(:);cpl(:);hph(:);cph(:);hsl(:);csl(:);hsh(:);csh(:)];
            if isempty(allv), allv = [0;1]; end
            mn = min(allv); mx = max(allv);
            pad = 0.05*(mx-mn + eps);
            ylim([mn-pad, mx+pad]);
        else
            ylim(yLimitLocal);
        end

        yl = get(gca,'YLim');

        % Labels and divider
        text(1.5, yl(2)*0.92, 'Pressure', 'FontSize', 10, 'HorizontalAlignment','center', 'Color','k');
        text(4.5, yl(2)*0.92, 'Sound',    'FontSize', 10, 'HorizontalAlignment','center', 'Color','k');
        line([3 3], yl, 'Color',[0.8 0.8 0.8], 'LineStyle','--', 'LineWidth', 1.5);

        xticks([1 2 4 5]);
        xticklabels({'Low','High','Low','High'});
        set(gca, 'FontSize', 10);
    end
end
