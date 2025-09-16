function Rall = lmm_baseline_min(T)
% Minimal LMMs: HC vs CBP (Low/High), random intercept per subject.

% Types
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
        lme = fitlme(Tk,'value ~ GroupBin*intensity + (1|subID)','FitMethod','REML');
        AT  = anova(lme,'DFMethod','Satterthwaite');

        % extract rows (use NaN when missing)
        rowG  = strcmp(AT.Term,'GroupBin');
        rowGI = strcmp(AT.Term,'GroupBin:intensity');

        FG   = iff(any(rowG),  AT.FStat(rowG),  NaN);
        df1G = iff(any(rowG),  AT.DF1(rowG),    NaN);
        df2G = iff(any(rowG),  AT.DF2(rowG),    NaN);
        pG   = iff(any(rowG),  AT.pValue(rowG), NaN);

        FGI  = iff(any(rowGI), AT.FStat(rowGI), NaN);
        pGI  = iff(any(rowGI), AT.pValue(rowGI),NaN);

        R = [R; {string(mk), height(Tk), FG, df1G, df2G, pG, FGI, pGI}];
    end

    % FDR on numeric, non-NaN entries
    if ~isempty(R)
        pg  = R.p_group;   pg  = double(pg);   maskG  = ~isnan(pg);
        pgi = R.p_gxint;   pgi = double(pgi);  maskGI = ~isnan(pgi);

        qG  = nan(size(pg));
        qGI = nan(size(pgi));

        if exist('fdr_bh','file')==2 && any(maskG)
            [~,~,~,qG(maskG)] = fdr_bh(pg(maskG));
        end
        if exist('fdr_bh','file')==2 && any(maskGI)
            [~,~,~,qGI(maskGI)] = fdr_bh(pgi(maskGI));
        end

        R.q_group  = qG;
        R.q_gxint  = qGI;
    end

    Rall.(mod) = R;
end
end

function out = iff(cond, a, b)
% inline ternary returning scalars (avoids [] creeping in)
if cond, out = a; else, out = b; end
end
