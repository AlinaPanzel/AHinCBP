function neural_baseline = get_neural_baseline(d, lo)
% Baseline-only (S1) neural table incl. HC and CBP, with ROIs + MVPA.
% Returns long format:
%   subID | timepoint | group | intensity | modality | measure | value
%
% Uses:
%   lo.roi.<G>.S1.<ROI>.(all_s_l/all_s_h/all_t_l/all_t_h)
%   lo.fm_mvpa.<G>.S1.(all_s_l/all_s_h/all_t_l/all_t_h).FM_PAIN / FM_MSS
%   lo.na_mvpa.<G>.S1.(all_s_l/all_s_h/all_t_l/all_t_h).general / sound
%
% where <G> is G1,G2,G3,HC.

%% Column headers: ROIs + MVPA (must match every row below)
headers = {'ID','Time','Group', ...
           'A1','MGN','IC','S1_L','S1_R', ...
           'm_ventral_insula','m_dorsal_insula','m_posterior_insula', ...
           'mPFC','precuneus','PCC', ...
           'FM_PAIN','FM_MSS','general','sound'};

% Helper to build one (group x modality x intensity) block
    function M = block(Gtag, grpCode, modeTag, intTag)
        % Gtag  : 'G1','G2','G3','HC'
        % grpCode: 1,2,3,0
        % modeTag: 's' (sound) or 't' (pressure/thumb)
        % intTag : 'l' (low) or 'h' (high)
        sess = 'S1';
        if modeTag=='s'
            fstr = sprintf('all_s_%c',intTag);
        else
            fstr = sprintf('all_t_%c',intTag);
        end

        % IDs/Time
        ids  = d.(Gtag).S1.id;
        time = d.(Gtag).S1.time;

        % ROI pulls
        A1  = lo.roi.(Gtag).(sess).A1.(fstr);
        MGN = lo.roi.(Gtag).(sess).MGN.(fstr);
        IC  = lo.roi.(Gtag).(sess).IC.(fstr);
        S1L = lo.roi.(Gtag).(sess).S1_L.(fstr);
        S1R = lo.roi.(Gtag).(sess).S1_R.(fstr);
        mVI = lo.roi.(Gtag).(sess).m_ventral_insula.(fstr);
        mDI = lo.roi.(Gtag).(sess).m_dorsal_insula.(fstr);
        mPI = lo.roi.(Gtag).(sess).m_posterior_insula.(fstr);
        mPFC= lo.roi.(Gtag).(sess).mPFC.(fstr);
        pre = lo.roi.(Gtag).(sess).precuneus.(fstr);
        PCC = lo.roi.(Gtag).(sess).PCC.(fstr);

        % MVPA pulls
        FMp  = lo.fm_mvpa.(Gtag).(sess).(fstr).FM_PAIN;
        FMss = lo.fm_mvpa.(Gtag).(sess).(fstr).FM_MSS;
        NAgen= lo.na_mvpa.(Gtag).(sess).(fstr).general;
        NAsnd= lo.na_mvpa.(Gtag).(sess).(fstr).sound;

        % Pack
        M = [ids, time, grpCode*ones(numel(ids),1), ...
             A1, MGN, IC, S1L, S1R, mVI, mDI, mPI, mPFC, pre, PCC, ...
             FMp, FMss, NAgen, NAsnd];
    end

%% ---------- Build per-modality × intensity wide tables ----------
% SOUND low/high
S_l = [block('G1',1,'s','l'); block('G2',2,'s','l'); block('G3',3,'s','l'); block('HC',0,'s','l')];
S_h = [block('G1',1,'s','h'); block('G2',2,'s','h'); block('G3',3,'s','h'); block('HC',0,'s','h')];

T_s_l = array2table(S_l, 'VariableNames', headers); T_s_l = sortrows(T_s_l,'ID');
T_s_h = array2table(S_h, 'VariableNames', headers); T_s_h = sortrows(T_s_h,'ID');

% PRESSURE low/high
T_l = [block('G1',1,'t','l'); block('G2',2,'t','l'); block('G3',3,'t','l'); block('HC',0,'t','l')];
T_h = [block('G1',1,'t','h'); block('G2',2,'t','h'); block('G3',3,'t','h'); block('HC',0,'t','h')];

T_t_l = array2table(T_l, 'VariableNames', headers); T_t_l = sortrows(T_t_l,'ID');
T_t_h = array2table(T_h, 'VariableNames', headers); T_t_h = sortrows(T_t_h,'ID');

%% ---------- Add Intensity + Modality keys ----------
addKM = @(T,intVal,modStr) movevars( ...
                 addvars(T, intVal*ones(height(T),1), repmat(string(modStr),height(T),1), ...
                         'NewVariableNames',{'Intensity','Modality'}), ...
                 {'Intensity','Modality'}, 'After','Group');

T_s_l = addKM(T_s_l,1,'Sound');    T_s_h = addKM(T_s_h,2,'Sound');
T_t_l = addKM(T_t_l,1,'Pressure'); T_t_h = addKM(T_t_h,2,'Pressure');

%% ---------- Align, concat, and melt to long ----------
commonVars = T_s_l.Properties.VariableNames;
T_s_h = T_s_h(:,commonVars);
T_t_l = T_t_l(:,commonVars);
T_t_h = T_t_h(:,commonVars);

W = [T_s_l; T_s_h; T_t_l; T_t_h];
W = sortrows(W, {'ID','Time'});

keyVars  = {'ID','Time','Group','Intensity','Modality'};
dataVars = setdiff(W.Properties.VariableNames, keyVars, 'stable');

neural_baseline = stack(W, dataVars, ...
                            'NewDataVariableName','value', ...
                            'IndexVariableName','measure');

% Final tidy types/labels
neural_baseline.Properties.VariableNames{'ID'}        = 'subID';
neural_baseline.Properties.VariableNames{'Time'}      = 'timepoint';
neural_baseline.Properties.VariableNames{'Group'}     = 'group';
neural_baseline.Properties.VariableNames{'Intensity'} = 'intensity';

neural_baseline.modality  = categorical(string(neural_baseline.Modality));
neural_baseline.Modality  = [];
neural_baseline.intensity = categorical(neural_baseline.intensity,[1 2],{'Low','High'});

neural_baseline = sortrows(neural_baseline, ...
    {'subID','timepoint','group','modality','intensity','measure'});

% After creating neural_baseline_roi
neural_baseline.GroupBin = categorical( ...
    ismember(neural_baseline.group, [1 2 3]), ...
    [0 1], {'HC','CBP'});

neural_baseline.timepoint = [];
% Rename first
neural_baseline.Properties.VariableNames{'group'}    = 'treatgroup';
neural_baseline.Properties.VariableNames{'GroupBin'} = 'group';

% Now put the columns in the exact order you want
neural_baseline = neural_baseline(:, ...
    {'subID','group','treatgroup','measure','modality','intensity','value'});

% (Optional) also sort the rows by these keys
% I'd usually omit 'value' here to avoid reordering by magnitude.
neural_baseline = sortrows(neural_baseline, ...
    {'subID','group','treatgroup','measure','modality','intensity'});

end
