function R = lmm_neural_baseline_allmeasures(neural_baseline)
% Baseline LMMs (S1) for all measures (ROIs + MVPA), SOUND & PRESSURE if present
% Compares HC (group=0) vs CBP (group in 1:3), with Intensity (1=Low,2=High)
% Fixed: GroupBin * Intensity; Random: (1|subID)
% Robust to empty cells / unused levels; per-group outlier removal

T = neural_baseline;
T = T(T.timepoint==1, :);                         % baseline only

% Types
T.value     = double(T.value);
T.subID     = categorical(T.subID);
T.intensity = categorical(T.intensity,[1 2],{'Low','High'});
T.measure   = categorical(T.measure);
T.group     = double(T.group);

% Binary group: HC vs CBP (any of 1–3)
T.GroupBin  = categorical(ismember(T.group, [1 2 3]), [0 1], {'HC','CBP'});

measures = categories(T.measure);

R = table( ...
    string.empty, zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    nan(0,1), nan(0,1), ...  % means
    string.empty, ...
    'VariableNames', {'measure','N','nHC','nCBP','nLow','nHigh', ...
                      'F_group','df1_group','df2_group','p_group', ...
                      'F_gxint','p_gxint', ...
                      'formula'});

fprintf('=== Baseline LMMs (HC vs CBP) across measures ===\n');

for k = 1:numel(measures)
    m  = measures{k};
    Tk = T(T.measure==m, :);
    if height(Tk) < 8, continue; end

    % ---- per-group outlier removal (rmoutliers) ----
    [gID, gLabs] = findgroups(Tk.GroupBin);
    drop = false(height(Tk),1);
    for gi = 1:numel(gLabs)
        idx = (gID==gi);
        if any(idx)
            [~, out] = rmoutliers(Tk.value(idx));   % by median/MAD
            drop(find(idx,1,'first')-1+(1:sum(idx))) = out; 
        end
    end
    Tk(drop,:) = [];
    if height(Tk) < 6, continue; end

    % Drop unused levels to avoid singular design
    Tk.GroupBin  = removecats(Tk.GroupBin);
    Tk.intensity = removecats(Tk.intensity);

    lvG = categories(Tk.GroupBin);
    lvI = categories(Tk.intensity);

    % Count cells to see if the 2x2 is complete
    Tk.ones = ones(height(Tk),1);
    C = groupsummary(Tk, {'GroupBin','intensity'}, 'sum', 'ones');  % then use C.sum_ones
    hasG = numel(lvG)==2;
    hasI = numel(lvI)==2;
    full2x2 = hasG && hasI && height(C)==4 && all(C.GroupCount>0);

    % Choose estimable formula
    if full2x2
        form = 'value ~ GroupBin*intensity + (1|subID)';
        % To use a random slope for intensity (only if both levels per subject):
        % form = 'value ~ GroupBin*intensity + (1 + intensity | subID)';
    elseif hasG && hasI
        form = 'value ~ GroupBin + intensity + (1|subID)';
    elseif hasG
        form = 'value ~ GroupBin + (1|subID)';
    elseif hasI
        form = 'value ~ intensity + (1|subID)';
    else
        continue;
    end

    % Fit
    lme = fitlme(Tk, form, 'FitMethod','REML');

    % Pull group and interaction (if present)
    cn = string(lme.CoefficientNames);

    [FG,df1G,df2G,pG] = deal(NaN);
    ixG = find(cn=="GroupBin_CBP",1);
    if ~isempty(ixG)
        Lg = zeros(1,numel(cn)); Lg(ixG)=1;
        [pG,FG,df1G,df2G] = coefTest(lme, Lg, 0, 'DFMethod','Satterthwaite');
    end

    [FGI,pGI] = deal(NaN);
    ixGI = find(cn=="GroupBin_CBP:intensity_High",1);
    if ~isempty(ixGI)
        Li = zeros(1,numel(cn)); Li(ixGI)=1;
        [pGI,FGI] = coefTest(lme, Li, 0, 'DFMethod','Satterthwaite');
    end

    % Simple counts & means (post-outlier)
    nHC   = sum(Tk.GroupBin=="HC");
    nCBP  = sum(Tk.GroupBin=="CBP");
    nLow  = sum(Tk.intensity=="Low");
    nHigh = sum(Tk.intensity=="High");

    R = [R; {string(m), height(Tk), nHC, nCBP, nLow, nHigh, ...
             FG, df1G, df2G, pG, ...
             FGI, pGI, ...
             string(form)}];
end

% BH–FDR across measures
if ~isempty(R)
    [~,~,~,qG]  = fdr_bh(R.p_group);
    [~,~,~,qGI] = fdr_bh(R.p_gxint);
    R.q_group  = qG;
    R.q_gxint  = qGI;

    % Nicely rounded print
    RR = R;
    RR.F_group = round(RR.F_group,3); RR.df1_group = round(RR.df1_group,2); RR.df2_group = round(RR.df2_group,2);
    RR.p_group = round(RR.p_group,4); RR.q_group   = round(RR.q_group,4);
    RR.F_gxint = round(RR.F_gxint,3); RR.p_gxint   = round(RR.p_gxint,4); RR.q_gxint = round(RR.q_gxint,4);

    RR = sortrows(RR, {'p_group','p_gxint'});
    disp('=== LMM (baseline) per measure: Group main effect and Group×Intensity ===');
    disp(RR(:, {'measure','N','nHC','nCBP','F_group','df1_group','df2_group','p_group','q_group','F_gxint','p_gxint','q_gxint','formula'}));
else
    disp('No estimable models.');
end
end
