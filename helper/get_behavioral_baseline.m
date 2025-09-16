function data_baseline = get_behavioral_baseline(d) 
%GET_BASELINE_TABLE Build a long table for baseline (S1) behavioral ratings.
% Returns one table with columns:
%   ID, Gender, Age, Group (Healthy/CBP), Modality (Pressure/Sound),
%   Intensity (Low/High), Rating
%
% Assumes d.metadata has (at minimum):
%   id, time, gender, age,
%   acute_mean_thumb_lo, acute_mean_thumb_hi,
%   acute_mean_sound_lo, acute_mean_sound_hi,
%   (either) is_patient OR group (with CBP in [1 2 3], HC otherwise).

    % ---- baseline subset ----
    M = d.metadata;
    M = M(M.time == 1, :);              % baseline only

    % ---- group coding → 'Healthy' / 'CBP' ----
    if ismember('is_patient', M.Properties.VariableNames)
        isCBP = logical(M.is_patient);
    elseif ismember('group', M.Properties.VariableNames)
        isCBP = ismember(M.group, [1 2 3]);  % adjust if your codes differ
    else
        error('get_baseline_table: Cannot infer group. Need is_patient or group.');
    end
    Group = categorical(isCBP, [0 1], {'Healthy','CBP'});

    % ---- pull core columns (safe) ----
    ID     = M.id;
    Gender = M.gender;     % assumed 1=male, 2=female
    Age    = M.age;

    % ---- ratings (behavioral means) ----
    req = {'acute_mean_thumb_lo','acute_mean_thumb_hi', ...
           'acute_mean_sound_lo','acute_mean_sound_hi'};
    missing = setdiff(req, M.Properties.VariableNames);
    if ~isempty(missing)
        error('get_baseline_table: Missing columns in d.metadata: %s', strjoin(missing, ', '));
    end

    % Make sure these are column vectors
    pp_lo = M.acute_mean_thumb_lo(:);
    pp_hi = M.acute_mean_thumb_hi(:);
    aa_lo = M.acute_mean_sound_lo(:);
    aa_hi = M.acute_mean_sound_hi(:);

    % ---- build two wide tables (Pressure / Sound) with shared IDs ----
    T_pressure = table(ID, Gender, Age, Group, pp_lo, pp_hi, ...
                       'VariableNames', {'ID','Gender','Age','Group','Avg_lo','Avg_hi'});

    T_sound    = table(ID, Gender, Age, Group, aa_lo, aa_hi, ...
                       'VariableNames', {'ID','Gender','Age','Group','Avg_lo','Avg_hi'});

    % ---- reshape to long ----
    P_long = stack(T_pressure, {'Avg_lo','Avg_hi'}, ...
                   'NewDataVariableName','Rating', ...
                   'IndexVariableName','Intensity');

    S_long = stack(T_sound, {'Avg_lo','Avg_hi'}, ...
                   'NewDataVariableName','Rating', ...
                   'IndexVariableName','Intensity');

    % Label modalities
    P_long.Modality = categorical(repmat("Pressure", height(P_long), 1));
    S_long.Modality = categorical(repmat("Sound",    height(S_long), 1));

    % Clean up Intensity labels
    P_long.Intensity = categorical(strrep(string(P_long.Intensity), 'Avg_', ''), ...
                                   {'lo','hi'}, {'Low','High'});
    S_long.Intensity = categorical(strrep(string(S_long.Intensity), 'Avg_', ''), ...
                                   {'lo','hi'}, {'Low','High'});

    % Concatenate
    data_baseline = [P_long; S_long];

    % Optional: map gender to categorical
    if isnumeric(Gender) || islogical(Gender)
        % assumes 1=male, 2=female — adjust if different
        data_baseline.nGender = categorical(data_baseline.Gender, [1 2], {'male','female'});
    else
        data_baseline.nGender = categorical(string(data_baseline.Gender));
    end

    % Centered age (handy for LMM)
    data_baseline.cAge = data_baseline.Age - nanmean(data_baseline.Age);

    % Housekeeping: order columns nicely
    data_baseline = movevars(data_baseline, ...
        {'Modality','Intensity','Rating','Group','nGender','cAge'}, 'After','Age');

    % ---- sort by ID ----
    data_baseline = sortrows(data_baseline, "ID");
end
