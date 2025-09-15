function results = lmm_cbp_vs_hc_baseline(T)
% CBP (groups 1–3 pooled) vs HC (group 0) at BASELINE, per ROI and modality.
% Requires columns: subID, timepoint (1), group, intensity, measure, value, modality

% --- baseline only ---
if ismember('timepoint',T.Properties.VariableNames)
    T = T(T.timepoint==1,:);
elseif ismember('Time',T.Properties.VariableNames)
    T = T(T.Time==1,:);
else
    error('Need timepoint (or Time) with baseline==1.');
end

% --- normalize types/labels ---
need = {'subID','group','intensity','measure','value','modality'};
miss = setdiff(need, T.Properties.VariableNames);
if ~isempty(miss), error('Missing columns: %s', strjoin(miss,', ')); end

T.subID     = categorical(T.subID);
T.intensity = categorical(T.intensity);          % 1/2
T.measure   = string(T.measure);
T.modality  = categorical(lower(string(T.modality))); % 'pressure'/'sound'
T = T(~isnan(T.value) & isfinite(T.value), :);

% pool CBP (groups 1–3) vs HC (0)
isHC  = (T.group==0) | strcmpi(string(T.group),'hc');
T.isCBP = double(~isHC);

mods     = unique(T.modality);              % modalities present
measures = unique(string(T.measure));       % ROI list (strings)

results = table(string([]), string([]), ...
                nan(0,1),nan(0,1),nan(0,1),nan(0,1), ...
                nan(0,1),nan(0,1),nan(0,1),nan(0,1), ...
                'VariableNames',{'measure','modality', ...
                                 'F_group','df1_group','df2_group','p_group', ...
                                 'F_gxint','df1_gxint','df2_gxint','p_gxint'});

for m = 1:numel(mods)
    M  = mods(m);
    TM = T(T.modality==M,:);

    for k = 1:numel(measures)
        R  = measures(k);
        % -------- FIXED LINE: filter with TM.measure (not T.measure) --------
        Tk = TM(string(TM.measure) == R, :);
        % -------------------------------------------------------------------

        if height(Tk) < 6, continue; end

        try
            lme = fitlme(Tk, 'value ~ isCBP * intensity + (1|subID)', 'FitMethod','REML');
        catch ME
            warning('LMM failed (ROI=%s, modality=%s): %s', R, string(M), ME.message);
            continue
        end

        cn = string(lme.CoefficientNames);

        % main effect: isCBP
        iG = find(cn=="isCBP",1);
        if ~isempty(iG), L=zeros(1,numel(cn)); L(iG)=1; [pG,FG,df1G,df2G]=coefTest(lme,L);
        else,            [FG,df1G,df2G,pG]=deal(NaN);
        end

        % interaction: isCBP:intensity*
        J = find(contains(cn,"isCBP:intensity"));
        if ~isempty(J)
            L = zeros(numel(J), numel(cn));
            for j=1:numel(J), L(j,J(j))=1; end
            [pGI,FGI,df1GI,df2GI] = coefTest(lme,L);
        else
            [FGI,df1GI,df2GI,pGI]=deal(NaN);
        end

        results = [results; {R, string(M), FG, df1G, df2G, pG, FGI, df1GI, df2GI, pGI}];
    end
end

% BH–FDR within each modality for p_group
results.q_group = nan(height(results),1);
for m = 1:numel(mods)
    rows = results.modality==string(mods(m));
    p = results.p_group(rows); keep = ~isnan(p);
    if any(keep), [~,~,~,q] = fdr_bh(p(keep)); results.q_group(rows) = NaN; results.q_group(rows(keep)) = q; end
end

% tidy print
if ~isempty(results)
    Rshow = results;
    Rshow.F_group   = round(Rshow.F_group,3);
    Rshow.df1_group = round(Rshow.df1_group,2);
    Rshow.df2_group = round(Rshow.df2_group,2);
    Rshow.p_group   = round(Rshow.p_group,4);
    Rshow.q_group   = round(Rshow.q_group,4);
    disp('=== CBP vs HC at baseline (per ROI, per modality) ===');
    disp(Rshow(:,{'measure','modality','F_group','df1_group','df2_group','p_group','q_group'}));
else
    fprintf('[diag] No analyzable rows after filtering. Mods present: %s\n', strjoin(unique(string(T.modality))',', '));
end
end

% ---- BH helper ----
function [h, crit_p, adj_ci_cvrg, q] = fdr_bh(pvals,q,method,report)
if nargin<2||isempty(q), q=.05; end
if nargin<3||isempty(method), method='pdep'; end
if nargin<4, report=''; end
p = pvals(:); [ps, idx] = sort(p);
V = numel(ps); I = (1:V)';
cV = strcmpi(method,'pdep')*1 + strcmpi(method,'dep')*sum(1./(1:V)); if cV==0, cV=1; end
th = (I/V)*q/cV; isSig = ps<=th; cut = find(isSig,1,'last');
h = false(size(p)); if ~isempty(cut), h(idx(1:cut)) = true; end
crit_p = []; if ~isempty(cut), crit_p = ps(cut); end
q = nan(size(p)); if ~isempty(ps), q(idx) = min(1, cummin((V./I).*ps, 'reverse')); end
adj_ci_cvrg = NaN;
end
