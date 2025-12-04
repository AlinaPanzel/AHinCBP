function predict_treatment_changes(behavioral_longitudinal,d)

T = behavioral_longitudinal;   % auditory unpleasantness data

% keep only subjects with both pre and post in T
subList = unique(T.SubID);
keepT   = false(size(subList));

for i = 1:numel(subList)
    thisID = subList(i);
    tp = T.timepoint(T.SubID == thisID);
    keepT(i) = any(tp == 1) && any(tp == 2);
end

validSubs_T = subList(keepT);
T = T(ismember(T.SubID, validSubs_T), :);

M = d.metadata;   % clinical pain data

% keep only subjects with both pre and post in M
idList = unique(M.id);
keepM  = false(size(idList));

for i = 1:numel(idList)
    thisID = idList(i);
    tp = M.time(M.id == thisID);
    keepM(i) = any(tp == 1) && any(tp == 2);
end

validSubs_M = idList(keepM);
M = M(ismember(M.id, validSubs_M), :);

% keep only subjects that appear in both datasets
commonSubs = intersect(validSubs_T, validSubs_M);
T = T(ismember(T.SubID, commonSubs), :);
M = M(ismember(M.id,    commonSubs), :);

% baseline MSS from pre session, split by intensity
T_pre = T(T.timepoint == 1, :);

preMSS_byInt = varfun(@mean, T_pre, ...
    'InputVariables','rating', ...
    'GroupingVariables',{'SubID','intensity'});

MSS_wide = unstack(preMSS_byInt, 'mean_rating', 'intensity', ...
                   'VariableNamingRule','preserve');

vn = MSS_wide.Properties.VariableNames;
lowCol  = find(contains(vn,'1'));
highCol = find(contains(vn,'2'));

if isempty(lowCol) || isempty(highCol)
    error('Low or High intensity MSS column missing after unstack.');
end

MSS_wide.Properties.VariableNames{lowCol}  = 'MSS_low_pre';
MSS_wide.Properties.VariableNames{highCol} = 'MSS_high_pre';

subInfo = unique(T_pre(:, {'SubID','group'}));
MSS_tbl = join(subInfo, MSS_wide, 'Keys','SubID');

% pre and post clinical pain
M_pre  = M(M.time == 1, :);
M_post = M(M.time == 2, :);

M_pre  = renamevars(M_pre,  {'id','pain_avg'}, {'SubID','pain_pre'});
M_post = renamevars(M_post, {'id','pain_avg'}, {'SubID','pain_post'});

painTbl = join(M_pre(:,{'SubID','pain_pre'}), ...
               M_post(:,{'SubID','pain_post'}), ...
               'Keys','SubID');

% combine MSS and pain
predictTbl = join(MSS_tbl, painTbl, 'Keys','SubID');

% PRT indicator
predictTbl.isPRT = double(predictTbl.group == 1);

% post-tx pain = B1*pre-tx pain + B2*MSS_low + B3*(isPRT * MSS_low)

% low-intensity MSS model
disp('=== Low-intensity MSS model (pain_post) ===');
mdl_low = fitlm(predictTbl, ...
    'pain_post ~ pain_pre + MSS_low_pre + MSS_low_pre:isPRT');
disp(mdl_low);

% high-intensity MSS models
disp('=== High-intensity MSS model (pain_post) ===');
mdl_high = fitlm(predictTbl, ...
    'pain_post ~ pain_pre + MSS_high_pre + MSS_high_pre:isPRT');
disp(mdl_high);


end


