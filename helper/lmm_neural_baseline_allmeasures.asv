function Rall = lmm_neural_baseline_allmeasures(neural_baseline)
% Baseline LMMs (S1) for all measures (ROIs + MVPA), run separately for
% SOUND and PRESSURE modalities.
%
% Compares HC (group=0) vs CBP (group in 1:3), with Intensity (Low vs High).
% Fixed: GroupBin * Intensity; Random: (1|subID).
% Outliers removed within group before fitting.

T = neural_baseline_clean;
T = T(T.timepoint==1, :);   % baseline only

% Types
T.value = double(T.value);
T.subID = categorical(T.subID);

% --- FIX: make 'intensity' categorical robustly ---
if isnumeric(T.intensity)
    % already numeric codes 1/2
    T.intensity = categorical(T.intensity, [1 2], {'Low','High'});
else
    % text or mixed — normalize to "Low"/"High" then categorize
    si = lower(strtrim(string(T.intensity)));
    si(ismember(si, ["1","low"]))  = "Low";
    si(ismember(si, ["2","high"])) = "High";
    T.intensity = categorical(si, ["Low","High"]);  % valueset matches string class
end
% ---------------------------------------------------

T.measure  = categorical(T.measure);
T.group    = double(T.group);
T.modality = categorical(T.modality);

% Binary group: HC vs CBP
T.GroupBin = categorical(ismember(T.group,[1 2 3]), [0 1], {'HC','CBP'});

modalities = categories(T.modality);
Rall = struct();

for mIdx = 1:numel(modalities)
    mod = modalities{mIdx};
    fprintf('\n=== Baseline LMMs (HC vs CBP) for %s ===\n', upper(mod));

    Tm = T(T.modality==mod,:);
    measures = categories(Tm.measure);

    R = table( ...
        string.empty, zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
        nan(0,1), nan(0,1), ...
        string.empty, ...
        'VariableNames', {'measure','N','nHC','nCBP','nLow','nHigh', ...
        'F_group','df1_group','df2_group','p_group', ...
        'F_gxint','p_gxint', ...
        'formula'});

    for k = 1:numel(measures)
        meas = measures{k};
        Tk = Tm(Tm.measure==meas,:);
        if height(Tk) < 8, continue; end

        % Drop unused levels
        Tk.GroupBin  = removecats(Tk.GroupBin);
        Tk.intensity = removecats(Tk.intensity);

        lvG = categories(Tk.GroupBin);
        lvI = categories(Tk.intensity);

        % Count cells
        Tk.ones = ones(height(Tk),1);
        C = groupsummary(Tk, {'GroupBin','intensity'}, 'sum', 'ones');
        hasG = numel(lvG)==2;
        hasI = numel(lvI)==2;
        full2x2 = hasG && hasI && height(C)==4 && all(C.GroupCount>0);

        % Choose formula
        if full2x2
            form = 'value ~ GroupBin*intensity + (1|subID)';
        elseif hasG && hasI
            form = 'value ~ GroupBin + intensity + (1|subID)';
        elseif hasG
            form = 'value ~ GroupBin + (1|subID)';
        elseif hasI
            form = 'value ~ intensity + (1|subID)';
        else
            continue;
        end

        % Fit model
        lme = fitlme(Tk, form, 'FitMethod','REML');
        cn = string(lme.CoefficientNames);

        % Group effect
        [FG,df1G,df2G,pG] = deal(NaN);
        ixG = find(cn=="GroupBin_CBP",1);
        if ~isempty(ixG)
            Lg = zeros(1,numel(cn)); Lg(ixG)=1;
            [pG,FG,df1G,df2G] = coefTest(lme, Lg, 0, 'DFMethod','Satterthwaite');
        end

        % Interaction
        [FGI,pGI] = deal(NaN);
        ixGI = find(cn=="GroupBin_CBP:intensity_High",1);
        if ~isempty(ixGI)
            Li = zeros(1,numel(cn)); Li(ixGI)=1;
            [pGI,FGI] = coefTest(lme, Li, 0, 'DFMethod','Satterthwaite');
        end

        AT  = anova(lme,'DFMethod','Satterthwaite');      % table with Term, DF1, DF2, FStat, pValue
        rowG  = strcmp(AT.Term,'GroupBin');
        rowGI = strcmp(AT.Term,'GroupBin:intensity');

        [FG,df1G,df2G,pG] = deal(NaN);
        if any(rowG)
            FG   = AT.FStat(rowG);
            df1G = AT.DF1(rowG);
            df2G = AT.DF2(rowG);
            pG   = AT.pValue(rowG);
        end

        [FGI,pGI] = deal(NaN);
        if any(rowGI)
            FGI = AT.FStat(rowGI);
            pGI = AT.pValue(rowGI);
        end

        % Counts
        nHC   = sum(Tk.GroupBin=="HC");
        nCBP  = sum(Tk.GroupBin=="CBP");
        nLow  = sum(Tk.intensity=="Low");
        nHigh = sum(Tk.intensity=="High");

        % Store
        R = [R; {string(meas), height(Tk), nHC, nCBP, nLow, nHigh, ...
            FG, df1G, df2G, pG, FGI, pGI, string(form)}];
    end

    % FDR corrections
    if ~isempty(R)
        [~,~,~,qG]  = fdr_bh(R.p_group);
        [~,~,~,qGI] = fdr_bh(R.p_gxint);
        R.q_group  = qG;
        R.q_gxint  = qGI;

        % Round + print
        RR = R;
        RR.F_group = round(RR.F_group,3); RR.df1_group = round(RR.df1_group,2); RR.df2_group = round(RR.df2_group,2);
        RR.p_group = round(RR.p_group,4); RR.q_group   = round(RR.q_group,4);
        RR.F_gxint = round(RR.F_gxint,3); RR.p_gxint   = round(RR.p_gxint,4); RR.q_gxint = round(RR.q_gxint,4);

        RR = sortrows(RR, {'p_group','p_gxint'});
        disp(RR(:, {'measure','N','nHC','nCBP','F_group','df1_group','df2_group','p_group','q_group','F_gxint','p_gxint','q_gxint','formula'}));
    else
        disp('No estimable models.');
    end

    Rall.(mod) = R;
end
end
