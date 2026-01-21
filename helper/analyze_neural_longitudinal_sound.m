function [results, models] = analyze_neural_longitudinal_treatment(neural_longitudinal_sound)
% Fit & test treatment effects with LMM (random slopes & intercepts)
% Model: value ~ group*timepoint + intensity + (1 + timepoint | subID)
% Reports: PRT vs Placebo and PRT vs UC contrasts with b, 95% CI, p-values

T = neural_longitudinal_sound;

% ---- Set up & clean ----
T.value     = double(T.value);
T.subID     = categorical(T.subID);
T.timepoint = double(T.timepoint);   % 1=baseline, 2=post
T.intensity = double(T.intensity);   % 1=low, 2=high
T.measure   = categorical(string(T.measure));

% Set PRT as reference group for contrasts (PRT vs Placebo, PRT vs UC)
T.group = categorical(T.group, [1 2 3], {'PRT', 'Placebo', 'UC'});
T.group = reordercats(T.group, {'PRT', 'Placebo', 'UC'});

measures = categories(T.measure);

% ---- Preallocate results table ----
results = table('Size', [0, 12], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'NumObs', ...
                      'b_PRTvsPlacebo', 'CI_lo_Placebo', 'CI_hi_Placebo', 'p_PRTvsPlacebo', ...
                      'b_PRTvsUC', 'CI_lo_UC', 'CI_hi_UC', 'p_PRTvsUC', ...
                      'F_joint', 'p_joint'});

models = struct('measure', {}, 'lme', {});

fprintf('\n=== Neural Longitudinal Treatment Effects ===\n');
fprintf('Model: value ~ group*timepoint + intensity + (1 + timepoint | subID)\n');
fprintf('Reference group: PRT\n\n');

for k = 1:numel(measures)
    m  = measures{k};
    Tk = T(T.measure == m, :);
    
    % ---- Fit LMM: Simple model (pooled across intensity) ----
    lme = fitlme(Tk, 'value ~ group*timepoint + intensity + (1 + timepoint | subID)', ...
                 'FitMethod', 'REML');
    
    % ---- Extract coefficients ----
    coefTbl   = lme.Coefficients;
    coefNames = string(coefTbl.Name);
    
    % Find group × timepoint interaction terms
    idx_Placebo = find(contains(coefNames, 'Placebo') & contains(coefNames, 'timepoint'));
    idx_UC      = find(contains(coefNames, 'UC') & contains(coefNames, 'timepoint'));
    
    % PRT vs Placebo × Time
    if ~isempty(idx_Placebo)
        idx_Placebo = idx_Placebo(1);
        b_Placebo     = coefTbl.Estimate(idx_Placebo);
        se_Placebo    = coefTbl.SE(idx_Placebo);
        p_Placebo     = coefTbl.pValue(idx_Placebo);
        ci_lo_Placebo = b_Placebo - 1.96 * se_Placebo;
        ci_hi_Placebo = b_Placebo + 1.96 * se_Placebo;
    else
        b_Placebo = NaN; ci_lo_Placebo = NaN; ci_hi_Placebo = NaN; p_Placebo = NaN;
    end
    
    % PRT vs UC × Time
    if ~isempty(idx_UC)
        idx_UC = idx_UC(1);
        b_UC     = coefTbl.Estimate(idx_UC);
        se_UC    = coefTbl.SE(idx_UC);
        p_UC     = coefTbl.pValue(idx_UC);
        ci_lo_UC = b_UC - 1.96 * se_UC;
        ci_hi_UC = b_UC + 1.96 * se_UC;
    else
        b_UC = NaN; ci_lo_UC = NaN; ci_hi_UC = NaN; p_UC = NaN;
    end
    
    % ---- Joint F-test of group × timepoint ----
    J = find(contains(coefNames, "group") & contains(coefNames, "timepoint"));
    if ~isempty(J)
        L = zeros(numel(J), numel(coefNames));
        for i = 1:numel(J), L(i, J(i)) = 1; end
        [p_joint, F_joint, ~, ~] = coefTest(lme, L);
    else
        F_joint = NaN; p_joint = NaN;
    end
    
    % ---- Store results ----
    results = [results; {string(m), height(Tk), ...
                         b_Placebo, ci_lo_Placebo, ci_hi_Placebo, p_Placebo, ...
                         b_UC, ci_lo_UC, ci_hi_UC, p_UC, ...
                         F_joint, p_joint}];
    
    models(end+1).measure = string(m);
    models(end).lme = lme;
end

% ---- Print results ----
fprintf('%-20s | %-40s | %-40s\n', 'Measure', 'PRT vs Placebo × Time', 'PRT vs UC × Time');
fprintf('%s\n', repmat('-', 1, 105));

for i = 1:height(results)
    % Format PRT vs Placebo
    if results.p_PRTvsPlacebo(i) < 0.05
        str_plac = sprintf('b=%6.2f [%5.2f,%5.2f] p=%.3f *', ...
            results.b_PRTvsPlacebo(i), results.CI_lo_Placebo(i), results.CI_hi_Placebo(i), results.p_PRTvsPlacebo(i));
    else
        str_plac = sprintf('b=%6.2f [%5.2f,%5.2f] p=%.3f', ...
            results.b_PRTvsPlacebo(i), results.CI_lo_Placebo(i), results.CI_hi_Placebo(i), results.p_PRTvsPlacebo(i));
    end
    
    % Format PRT vs UC
    if results.p_PRTvsUC(i) < 0.05
        str_uc = sprintf('b=%6.2f [%5.2f,%5.2f] p=%.3f *', ...
            results.b_PRTvsUC(i), results.CI_lo_UC(i), results.CI_hi_UC(i), results.p_PRTvsUC(i));
    else
        str_uc = sprintf('b=%6.2f [%5.2f,%5.2f] p=%.3f', ...
            results.b_PRTvsUC(i), results.CI_lo_UC(i), results.CI_hi_UC(i), results.p_PRTvsUC(i));
    end
    
    fprintf('%-20s | %-40s | %-40s\n', results.Measure(i), str_plac, str_uc);
end

fprintf('%s\n', repmat('-', 1, 105));
fprintf('* = p < 0.05\n\n');

% ---- Highlight significant findings ----
fprintf('=== Significant Treatment Effects (p < 0.05) ===\n\n');

sig_found = false;
for i = 1:height(results)
    if results.p_PRTvsPlacebo(i) < 0.05
        fprintf('%s: PRT vs Placebo × Time: b = %.2f, 95%% CI [%.2f, %.2f], p = %.3f\n', ...
            results.Measure(i), results.b_PRTvsPlacebo(i), ...
            results.CI_lo_Placebo(i), results.CI_hi_Placebo(i), results.p_PRTvsPlacebo(i));
        sig_found = true;
    end
    if results.p_PRTvsUC(i) < 0.05
        fprintf('%s: PRT vs UC × Time: b = %.2f, 95%% CI [%.2f, %.2f], p = %.3f\n', ...
            results.Measure(i), results.b_PRTvsUC(i), ...
            results.CI_lo_UC(i), results.CI_hi_UC(i), results.p_PRTvsUC(i));
        sig_found = true;
    end
end

if ~sig_found
    fprintf('None.\n');
end

fprintf('\n');

end