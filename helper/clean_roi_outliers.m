function [neural_baseline_roi_clean, outlier_summary] = clean_roi_outliers(neural_baseline_roi)
% Clean ROI table by removing outliers within each ROI × GroupBin × Modality × Intensity.
% Returns the cleaned table and a summary of removals per cell.

T = neural_baseline_roi;

% ----- GroupBin (HC vs CBP) -----
T.GroupBin = categorical(ismember(T.group,[1 2 3]), [0 1], {'HC','CBP'});

% ----- Normalize INTENSITY to {'Low','High'} -----
% Accept numbers 1/2 or strings like 'low'/'high'
ival = string(T.intensity);
% If it looks numeric, convert:
try
    inum = double(T.intensity);
    if ~any(isnan(inum)), ival = string(inum); end
catch
    % ok as string already
end
isHigh = (ival=="2") | (lower(ival)=="high") | (lower(ival)=="hi") | (upper(ival)=="H");
isLow  = (ival=="1") | (lower(ival)=="low")  | (lower(ival)=="lo") | (upper(ival)=="L");
ival(~(isLow|isHigh)) = missing;  % anything else → missing
ival(isLow)  = "Low";
ival(isHigh) = "High";
T.intensity  = categorical(ival, ["Low","High"]);

% ----- Normalize MODALITY to {'Sound','Pressure'} -----
mval = string(T.modality);
isSound    = lower(mval)=="sound"    | lower(mval)=="audio";
isPressure = lower(mval)=="pressure" | lower(mval)=="thumb" | lower(mval)=="pp";
mval(~(isSound|isPressure)) = missing;
mval(isSound)    = "Sound";
mval(isPressure) = "Pressure";
T.modality       = categorical(mval, ["Sound","Pressure"]);

% ----- Ensure MEASURE is categorical -----
T.measure = categorical(string(T.measure));

% ----- Flag outliers within each (ROI × GroupBin × Modality × Intensity) cell -----
[grp, roiL, grpL, modL, intL] = findgroups(T.measure, T.GroupBin, T.modality, T.intensity);
isOut = false(height(T),1);

for k = 1:max(grp)
    idx = (grp == k);
    n = nnz(idx);
    if n >= 3   % need a few points to define outliers
        [~, outk] = rmoutliers(T.value(idx));
        isOut(idx) = outk;
    end
end

% Attach flag and drop outliers
T.outlier = isOut;
neural_baseline_roi_clean = T(~T.outlier, :);

% ----- Summarize removals per cell -----
outlier_summary = table(roiL, grpL, modL, intL, ...
    splitapply(@numel, T.value, grp), ...
    splitapply(@sum,   double(T.outlier), grp), ...
    'VariableNames', {'ROI','GroupBin','Modality','Intensity','N_total','N_removed'});

outlier_summary = sortrows(outlier_summary, {'ROI','Modality','Intensity','GroupBin'});

% ----- Optional: quick printout -----
fprintf('\n=== Outlier removal summary (per ROI × Modality × Intensity × Group) ===\n');
fprintf('%-16s %-6s %-9s %-6s %8s %10s\n','ROI','Group','Modality','Int','N_total','N_removed');
fprintf('%s\n', repmat('-',1,64));
for i = 1:height(outlier_summary)
    fprintf('%-16s %-6s %-9s %-6s %8d %10d\n', ...
        char(outlier_summary.ROI(i)), ...
        char(outlier_summary.GroupBin(i)), ...
        char(outlier_summary.Modality(i)), ...
        char(outlier_summary.Intensity(i)), ...
        outlier_summary.N_total(i), ...
        outlier_summary.N_removed(i));
end
fprintf('----------------------------------------------------------------\n');
fprintf('Total removed: %d of %d (%.1f%%)\n', ...
    nansum(outlier_summary.N_removed), nansum(outlier_summary.N_total), ...
    100*nansum(outlier_summary.N_removed)/max(1,nansum(outlier_summary.N_total)));
end
