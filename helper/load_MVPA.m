function lo = load_MVPA(lo)

Treatment = {'HC','G1','G2','G3'};
sessions  = {'S1','S2'};
field     = {'all_t_l','all_t_h','all_s_l','all_s_h','completers_t_l','completers_t_h','completers_s_l','completers_s_h'};

for g = 1:numel(Treatment)
    grp = Treatment{g};
    if ~isfield(lo, grp)
        fprintf('MVPA: skip group %s (missing in lo)\n', grp);
        continue
    end

    for s = 1:numel(sessions)
        ses = sessions{s};

        % Example : skip HC at S2 
        if strcmp(grp,'HC') && strcmp(ses,'S2')
            fprintf('MVPA: skip %s %s (HC baseline-only)\n', grp, ses);
            continue
        end

        if ~isfield(lo.(grp), ses) || isempty(lo.(grp).(ses))
            fprintf('MVPA: skip %s %s (session missing/empty)\n', grp, ses);
            continue
        end

        % Only iterate over fields that actually exist for this group/session
        existing = fieldnames(lo.(grp).(ses));
        useFields = intersect(field, existing, 'stable');
        if isempty(useFields)
            fprintf('MVPA: skip %s %s (no target fields present)\n', grp, ses);
            continue
        end

        for i = 1:numel(useFields)
            fld = useFields{i};

            % Guard: empty or non-existent container
            if ~isfield(lo.(grp).(ses), fld) || isempty(lo.(grp).(ses).(fld))
                fprintf('MVPA: skip %s %s %s (field missing/empty)\n', grp, ses, fld);
                continue
            end

            % Apply patterns with try/catch so one failure doesn't stop the batch
            try
                lo.fm_mvpa.(grp).(ses).(fld) = apply_FM_patterns(lo.(grp).(ses).(fld), 'cosine_similarity');
            catch ME
                fprintf('MVPA: apply_FM_patterns failed for %s %s %s: %s\n', grp, ses, fld, ME.message);
                lo.fm_mvpa.(grp).(ses).(fld) = [];  % or NaN
            end

            try
                lo.na_mvpa.(grp).(ses).(fld) = apply_multiaversive_mpa2_patterns(lo.(grp).(ses).(fld));
            catch ME
                fprintf('MVPA: apply_multiaversive_mpa2_patterns failed for %s %s %s: %s\n', grp, ses, fld, ME.message);
                lo.na_mvpa.(grp).(ses).(fld) = [];  % or NaN
            end
        end
    end
end

end 