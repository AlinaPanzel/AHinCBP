function Rall = lmm_neural_baseline(neural_baseline_clean)
% LMMs: HC vs CBP (Low/High), random intercept per subject.
% Returns per-modality tables with APA-style report strings for Group and Group×Intensity.

% Types
T = neural_baseline_clean;
T.value     = double(T.value);
T.subID     = categorical(T.subID);
T.measure   = categorical(T.measure);
T.modality  = categorical(T.modality);
T.GroupBin  = categorical(T.GroupBin);
T.intensity = categorical(T.intensity);

mods = categories(T.modality);
Rall = struct();

for m = 1:numel(mods)
    mod  = mods{m};
    Tm   = T(T.modality==mod,:);
    meas = categories(Tm.measure);

    R = table(string.empty, zeros(0,1), ...
              nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
              nan(0,1), nan(0,1), ...
        'VariableNames', {'measure','N','F_group','df1_group','df2_group','p_group','F_gxint','p_gxint'});

    % Pre-allocate APA text fields
    R.apa_group  = strings(0,1);
    R.apa_gxint  = strings(0,1);

    for k = 1:numel(meas)
        mk = meas{k};
        Tk = Tm(Tm.measure==mk,:);
        if height(Tk) < 8, continue; end

        % require full 2x2 with non-empty cells
        C = groupsummary(Tk, {'GroupBin','intensity'});
        if numel(categories(Tk.GroupBin))~=2 || numel(categories(Tk.intensity))~=2 ...
                || height(C)~=4 || any(C.GroupCount==0)
            continue;
        end

        % fit
        TkLow = Tk(Tk.intensity=="Low",:);
        lme   = fitlme(TkLow, 'value ~ GroupBin + (1|subID)', 'FitMethod','REML');
       %lme = fitlme(Tk, 'value ~ GroupBin + intensity + (1 + intensity|subID)', 'FitMethod','REML');
        AT  = anova(lme,'DFMethod','Satterthwaite');

        % extract rows
        rowG  = strcmp(AT.Term,'GroupBin');
        rowGI = strcmp(AT.Term,'GroupBin:intensity');

        FG   = iff(any(rowG),  AT.FStat(rowG),  NaN);
        df1G = iff(any(rowG),  AT.DF1(rowG),    NaN);
        df2G = iff(any(rowG),  AT.DF2(rowG),    NaN);
        pG   = iff(any(rowG),  AT.pValue(rowG), NaN);

        FGI  = iff(any(rowGI), AT.FStat(rowGI), NaN);
        pGI  = iff(any(rowGI), AT.pValue(rowGI),NaN);

        % build APA strings
        apaG  = makeAPA('Group', FG, df1G, df2G, pG);
        apaGI = makeAPA('Group × Intensity', FGI, 1, NaN, pGI); % df1=1 for interaction

        R = [R; {string(mk), height(Tk), FG, df1G, df2G, pG, FGI, pGI, apaG, apaGI}];
    end

    % FDR
    if ~isempty(R)
        pg  = double(R.p_group);  maskG  = ~isnan(pg);
        pgi = double(R.p_gxint);  maskGI = ~isnan(pgi);

        qG  = nan(size(pg));
        qGI = nan(size(pgi));

        if exist('fdr_bh','file')==2 && any(maskG)
            [~,~,~,qG(maskG)]  = fdr_bh(pg(maskG));
        end
        if exist('fdr_bh','file')==2 && any(maskGI)
            [~,~,~,qGI(maskGI)] = fdr_bh(pgi(maskGI));
        end

        R.q_group = qG; 
        R.q_gxint = qGI;
    end

    Rall.(mod) = R;
end
end

function out = iff(cond,a,b)
if cond, out=a; else, out=b; end
end

function txt = makeAPA(label,F,df1,df2,p)
% Return APA-style string, e.g. "Group effect: F(1, 45.23) = 6.21, p = .016"
if isnan(F) || isnan(p)
    txt = label + " effect: not estimable";
    return
end
if p < .001
    pStr = "< .001";
else
    pStr = sprintf("= %.3f", p);
    % APA wants 3 decimals, strip leading zero
    pStr = erase(pStr,"0");
end
txt = sprintf("%s effect: F(%.0f, %.2f) = %.2f, p %s", ...
    label, df1, df2, F, pStr);
end