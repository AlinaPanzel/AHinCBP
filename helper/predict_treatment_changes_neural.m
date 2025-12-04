function predict_treatment_changes_neural(neural_longitudinal_sound, d)

N = neural_longitudinal_sound;   % neural data
M = d.metadata;                  % clinical pain

% keep only subjects with both pre and post in N
subList = unique(N.subID);
keepN   = false(size(subList));

for i = 1:numel(subList)
    thisID = subList(i);
    tp = N.timepoint(N.subID == thisID);
    keepN(i) = any(tp == 1) && any(tp == 2);
end

validSubs_N = subList(keepN);
N = N(ismember(N.subID, validSubs_N), :);

idList = unique(M.id); %same for clinpaim
keepM  = false(size(idList));

for i = 1:numel(idList)
    thisID = idList(i);
    tp = M.time(M.id == thisID);
    keepM(i) = any(tp == 1) && any(tp == 2);
end

validSubs_M = idList(keepM);
M = M(ismember(M.id, validSubs_M), :);

% keep only subjects present in both datasets
commonSubs = intersect(validSubs_N, validSubs_M);
N = N(ismember(N.subID, commonSubs), :);
M = M(ismember(M.id,    commonSubs), :);

% clinical pain pre/post
M_pre  = M(M.time == 1, :);
M_post = M(M.time == 2, :);

M_pre  = renamevars(M_pre,  {'id','pain_avg'}, {'subID','pain_pre'});
M_post = renamevars(M_post, {'id','pain_avg'}, {'subID','pain_post'});

painTbl = join(M_pre(:,{'subID','pain_pre'}), ...
               M_post(:,{'subID','pain_post'}), ...
               'Keys','subID');

% ROIs & MVPA patterns 
roiList = { ...
    'A1', ...
    'm_ventral_insula', ...
    'm_dorsal_insula', ...
    'm_posterior_insula', ...
    'mPFC', ...
    'precuneus', ...
    'FM_PAIN', ...
    'FM_MSS', ...
    'general', ...
    'sound' ...
};
for intensity = [1 2]

    if intensity == 1
        intLabel = 'low';
    else
        intLabel = 'high';
    end

    % baseline neural values for this intensity
    N_pre = N(N.timepoint == 1 & N.intensity == intensity, :);

    % make wide table: one column per measure (ROI/pattern)
    preWide = unstack(N_pre, 'value', 'measure', ...
                      'VariableNamingRule','preserve');

    % add group info, then clinical pain
    subInfo = unique(N_pre(:, {'subID','group'}));
    tmpTbl  = join(subInfo, preWide, 'Keys','subID');
    tmpTbl  = join(tmpTbl, painTbl, 'Keys','subID');

    % PRT indicator (1 = PRT, 0 = control/other)
    tmpTbl.isPRT = double(tmpTbl.group_subInfo == 1);

    for r = 1:numel(roiList)
        roi = roiList{r};

        if ~ismember(roi, tmpTbl.Properties.VariableNames)
            fprintf('Skipping %s (%s): not in table.\n', roi, intLabel);
            continue;
        end

        fprintf('=== %s (%s intensity) model (pain_post) ===\n', roi, intLabel);

        % Model:
        % pain_post = B1*pain_pre + B2*ROI + B3*(ROI * isPRT)
        f1 = sprintf('pain_post ~ pain_pre + %s + %s:isPRT', roi, roi);
        mdl = fitlm(tmpTbl, f1);
        disp(mdl);
    end
end

end
