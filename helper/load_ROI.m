function lo = load_ROI(a, lo)

% Load ROIs 
areas    = {'A1','MGN','IC','S1_L','S1_R','m_ventral_insula','m_dorsal_insula','m_posterior_insula','mPFC','precuneus','PCC'};
field    = {'all_t_l','all_t_h','all_s_l','all_s_h','completers_t_l','completers_t_h','completers_s_l','completers_s_h'};
Treatment = {'HC','G1','G2','G3'};
sessions  = {'S1','S2'};

for n = 1:numel(Treatment)
    grp = Treatment{n};
    if ~isfield(lo, grp), continue; end   % group not present → skip

    for s = 1:numel(sessions)
        ses = sessions{s};

        % HC are baseline-only
        if strcmp(grp,'HC') && strcmp(ses,'S2'), continue; end
        if ~isfield(lo.(grp), ses), continue; end

        % Only process fields that actually exist for this group/session
        existing   = fieldnames(lo.(grp).(ses));
        useFields  = intersect(field, existing, 'stable');
        if isempty(useFields), continue; end

        for u = 1:numel(areas)
            roiName = areas{u};
            if ~isfield(a, roiName), continue; end

            for i = 1:numel(useFields)
                fld = useFields{i};
                if ~isfield(lo.(grp).(ses), fld) || isempty(lo.(grp).(ses).(fld)), continue; end

                r = extract_roi_averages(lo.(grp).(ses).(fld), a.(roiName));
                lo.roi.(grp).(ses).(roiName).(fld) = cat(2, r.dat);
            end
        end
    end
end


end 