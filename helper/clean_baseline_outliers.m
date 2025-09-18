function [neural_baseline_clean, outlier_summary] = clean_baseline_outliers(neural_baseline)
% Clean ROI+MVPA table by removing outliers within each ROI × Group × Modality × Intensity.

T = neural_baseline;

% ---- Coerce types (light-touch, based on your table) ----
T.value     = double(T.value);
T.group     = categorical(string(T.group));                % 'HC' / 'CBP'
T.modality  = categorical(string(T.modality), ["Sound","Pressure"]);
T.intensity = categorical(string(T.intensity), ["Low","High"]);
T.measure   = categorical(string(T.measure));

% ---- Grouping: ROI × Group × Modality × Intensity ----
[G, roiL, grpL, modL, intL] = findgroups(T.measure, T.group, T.modality, T.intensity);

isOut = false(height(T),1);
K = max(G);

for k = 1:K
    idx = (G == k);
    if nnz(idx) >= 3                 % need ≥3 to define outliers robustly
        [~, outk] = rmoutliers(T.value(idx));   % default MAD rule
        isOut(idx) = outk;
    end
end

T.outlier = isOut;
neural_baseline_clean = T(~T.outlier, :);

% ---- Summary per cell ----
N_total   = splitapply(@numel, T.value, G);
N_removed = splitapply(@sum, double(T.outlier), G);
Pct_removed = 100 * (N_removed ./ max(N_total,1));

outlier_summary = table(roiL, grpL, modL, intL, N_total, N_removed, Pct_removed, ...
    'VariableNames', {'ROI','Group','Modality','Intensity','N_total','N_removed','Pct_removed'});

outlier_summary = sortrows(outlier_summary, {'ROI','Modality','Intensity','Group'});

% ---- Pretty print ----
fprintf('\n=== Outlier removal summary (ROI × Group × Modality × Intensity) ===\n');
fprintf('%-18s %-6s %-9s %-6s %8s %10s %8s\n','ROI','Group','Mod','Int','N_total','N_removed','Pct%');
fprintf('%s\n', repmat('-',1,70));
for i = 1:height(outlier_summary)
    fprintf('%-18s %-6s %-9s %-6s %8d %10d %7.1f\n', ...
        char(outlier_summary.ROI(i)), ...
        char(outlier_summary.Group(i)), ...
        char(outlier_summary.Modality(i)), ...
        char(outlier_summary.Intensity(i)), ...
        outlier_summary.N_total(i), ...
        outlier_summary.N_removed(i), ...
        outlier_summary.Pct_removed(i));
end
fprintf('----------------------------------------------------------------------\n');
fprintf('Total removed: %d of %d (%.1f%%)\n', ...
    nansum(N_removed), nansum(N_total), 100*nansum(N_removed)/max(1,nansum(N_total)));
end
