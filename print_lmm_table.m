
%% ================= helper =================
function [RES] = print_lmm_table(Tall, modality, allROIs)
fprintf('\n==== %s (LMM; group main effect per ROI) ====\n', upper(char(modality)));
T = Tall(Tall.modality==modality, :);

% table to collect results
RES = table(string.empty, zeros(0,1), zeros(0,1), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'ROI','nHC','nCBP','F_group','df1','df2','p_group','qBH'});

% space for p-values to do FDR at the end
p_all = nan(numel(allROIs),1);

% first pass: fit all ROIs, store rows (q filled later)
rows = strings(numel(allROIs),1);
keep  = false(numel(allROIs),1);

for i = 1:numel(allROIs)
    R = allROIs{i};
    Tk = T(T.measure==R, :);
    if isempty(Tk), continue; end

    % --- outlier removal within each GroupBin ---
    [g, gl] = findgroups(Tk.GroupBin);
    tf = false(height(Tk),1);
    for k = 1:numel(gl)
        idx = (g==k);
        [~, out] = rmoutliers(Tk.value(idx));
        tf(idx) = out;
    end
    Tk(tf,:) = [];

    if height(Tk) < 6, continue; end

    % counts (unique subIDs per group after cleaning)
    nHC  = numel(unique(Tk.subID(Tk.GroupBin=="HC")));
    nCBP = numel(unique(Tk.subID(Tk.GroupBin=="CBP")));

    % --- LMM (reference: HC & Low) ---
    lme = fitlme(Tk, 'value ~ GroupBin*intensity + (1|subID)', 'FitMethod','REML');

    % Group main effect: coefficient for GroupBin_CBP
    % With our ref levels, the columns are:
    % (Intercept)  GroupBin_CBP  intensity_High  GroupBin_CBP:intensity_High
    L = [0 1 0 0];
    [p,F,DF1,DF2] = coefTest(lme, L, 0, 'DFMethod','Satterthwaite');

    p_all(i) = p;
    keep(i)  = true;

    % store (q will be printed later)
    rows(i) = sprintf('%-16s %6d %6d %9.2f %7.0f %7.2f %10.4g %6s', ...
        R, nHC, nCBP, F, DF1, DF2, p, stars(p));

    % stash in RES too (q filled later)
    RES = [RES; {string(R), nHC, nCBP, F, DF1, DF2, p, NaN}]; %#ok<AGROW>
end

% FDR across ROIs (only those we kept)
p_kept = p_all(keep);
q = nan(size(p_all));
if any(keep)
    [~,~,~,q(keep)] = fdr_bh(p_kept);
end

% pretty header
fprintf('%-16s %6s %6s %9s %7s %7s %10s %6s %8s\n', ...
    'ROI','nHC','nCBP','F_group','df1','df2','p','sig','q(BH)');
fprintf('%s\n', repmat('-',1,86));

% print rows with q + q-stars
krow = 0;
for i = 1:numel(allROIs)
    if ~keep(i), continue; end
    krow = krow + 1;
    qtxt = sprintf('%.4g', q(i));
    fprintf('%s %8s%s\n', rows(i), qtxt, padstars(q(i)));
    RES.qBH(krow) = q(i);
end
end

function s = stars(val)
if isnan(val), s = ''; return; end
if val < 0.001, s = '***';
elseif val < 0.01, s = '**';
elseif val < 0.05, s = '*';
else, s = '';
end
end

function s = padstars(q)
% right-aligned star column for q
if isnan(q), s = ''; return; end
if q < 0.001, s = ' ***';
elseif q < 0.01, s = ' **';
elseif q < 0.05, s = ' *';
else, s = ' ';
end
end
