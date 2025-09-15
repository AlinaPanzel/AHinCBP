function behavioral_longitudinal = get_behavioral_longitudinal(d)
% Build long behavioral table: auditory unpleasantness ratings (CBP only).
% Returns table with:
%   SubID | timepoint(1,2) | group(1=PRT,2=Placebo,3=UC) |
%   intensity(1=Low,2=High) | rating | age | cAge | gender(1,2)

groups   = {'G1','G2','G3'};
sessions = {'S1','S2'};
grpCode  = struct('G1',1,'G2',2,'G3',3);  % numeric group codes

W = table();  % accumulator

for g = 1:numel(groups)
    grp = groups{g};
    for s = 1:numel(sessions)
        ses = sessions{s};
        if ~isfield(d, grp) || ~isfield(d.(grp), ses), continue; end
        S = d.(grp).(ses);
        if ~istable(S) || height(S)==0, continue; end

        % make sure all needed vars exist
        needVars = {'id','time','gender','age','acute_mean_sound_lo','acute_mean_sound_hi'};
        miss = setdiff(needVars, S.Properties.VariableNames);
        if ~isempty(miss)
            error('Missing columns in d.%s.%s: %s', grp, ses, strjoin(miss, ', '));
        end

        % build wide block from this session only
        T = table( ...
            S.id(:), ...
            S.time(:), ...
            repmat(grpCode.(grp), height(S), 1), ...
            S.gender(:), ...
            S.age(:), ...
            S.acute_mean_sound_lo(:), ...
            S.acute_mean_sound_hi(:), ...
            'VariableNames', {'ID','Time','group','gender','age','Avg_lo','Avg_hi'} ...
        );

        W = [W; T]; %#ok<AGROW>
    end
end

if isempty(W)
    behavioral_longitudinal = table();
    warning('get_behavioral_longitudinal: No rows found.');
    return
end

% center age
cAge = W.age - nanmean(W.age);

% reshape to long
L = stack(W, {'Avg_lo','Avg_hi'}, ...
          'NewDataVariableName','rating', ...
          'IndexVariableName','Intensity');

% numeric intensity
IntensityStr = string(L.Intensity);
intensity = nan(height(L),1);
intensity(IntensityStr=="Avg_lo") = 1;
intensity(IntensityStr=="Avg_hi") = 2;

% final table
behavioral_longitudinal = table( ...
    double(L.ID), ...
    double(L.Time), ...
    double(L.group), ...
    double(intensity), ...
    L.rating, ...
    L.age, ...
    double(L.gender), ...
    'VariableNames', {'SubID','timepoint','group','intensity','rating','age','gender'} ...
);

behavioral_longitudinal = sortrows(behavioral_longitudinal, {'SubID','timepoint','group','intensity'});

% Just to double check & print

% Count unique participants per group × timepoint
[G, groupID, timeID] = findgroups(behavioral_longitudinal.group, ...
                                  behavioral_longitudinal.timepoint);

nSubs = splitapply(@(x) numel(unique(x)), behavioral_longitudinal.SubID, G);

summaryTable = table(groupID, timeID, nSubs, ...
    'VariableNames', {'Group','Timepoint','N_Subjects'});
disp(summaryTable)

end
