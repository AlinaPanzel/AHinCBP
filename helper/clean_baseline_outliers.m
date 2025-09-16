function [neural_baseline_clean, outlier_summary] = clean_baseline_outliers(neural_baseline)
% Clean ROI+MVPA table by removing outliers within each ROI × GroupBin × Modality × Intensity.
% Returns the cleaned table and a summary of removals per cell.

T = neural_baseline;

% ----- GroupBin (HC vs CBP) -----
T.GroupBin = categorical(ismember(T.group,[1 2 3]), [0 1], {'HC','CBP'});

% ----- Normalize INTENSITY to {'Low','High'} -----
ival = string(T.intensity);
isHigh = (ival=="2") | strcmpi(ival,"high") | strcmpi(ival,"hi") | strcmpi(ival,"H");
isLow  = (ival=="1") | strcmpi(ival,"low")  | strcmpi(ival,"lo") | strcmpi(ival,"L");
ival(isLow)  = "Low";
ival(isHigh) = "High";
ival(~(isLow|isHigh)) = missing;
T.intensity  = categorical(ival, ["Low","High"]);

% ----- Normalize MODALITY to {'Sound','Pressure'} -----
mval = string(T.modality);
isSound    = strcmpi(mval,"sound")    | strcmpi(mval,"audio");
isPressure = strcmpi(mval,"pressure") | strcmpi(mval,"thumb") | strcmpi(mval,"pp");
mval(isSound)    = "Sound";
mval(isPressure) = "Pressure";
mval(~(isSound|isPressure)) = missing;
T.modality = categorical(mval, ["Sound","Pressure"]);

% ----- Ensure MEASURE is categorical (ROI + MVPA) -----
T.measure = categorical(string(T.measure));

% ----- Flag outliers within each (ROI × GroupBin × Modality × Intensity) -----
[G,roiL,grpL,modL,intL] = findgroups(T.measure, T.GroupBin, T.modality, T.intensity);
isOut = false(height(T),1);

for k = 1:max(G)
    idx = (G == k);
    if nnz(idx) >= 3   % need ≥3 values to define outliers
        [~, outk] = rmoutliers(T.value(idx));
        isOut(idx) = outk;
    end
end

T.outlier = isOut;
neural_baseline_clean = T(~T.outlier, :);

% ----- Summarize removals -----
N_total   = splitapply(@numel, T.value, G);
N_removed = splitapply(@sum, double(T.outlier), G);
perc      = 100 * (N_removed ./ max(N_total,1));

outlier_summary = table(roiL, grpL, modL, intL, N_total, N_removed, perc, ...
    'VariableNames', {'ROI','GroupBin','Modality','Intensity','N_total','N_removed','Pct_removed'});

outlier_summary = sortrows(outlier_summary, {'ROI','Modality','Intensity','GroupBin'});

% ----- Print summary -----
fprintf('\n=== Outlier removal summary (ROI × Group × Modality × Intensity) ===\n');
fprintf('%-18s %-6s %-9s %-6s %8s %10s %8s\n','ROI','Group','Mod','Int','N_total','N_removed','Pct%');
fprintf('%s\n', repmat('-',1,70));
for i = 1:height(outlier_summary)
    fprintf('%-18s %-6s %-9s %-6s %8d %10d %7.1f\n', ...
        char(outlier_summary.ROI(i)), ...
        char(outlier_summary.GroupBin(i)), ...
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
