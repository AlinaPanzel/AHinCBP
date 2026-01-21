function check_model_assumptions(mdl, T, mssVarName)

% mdl        : LinearModel from fitlm
% T          : table used to fit mdl (must contain isPRT and MSS variable)
% mssVarName : name of MSS variable in T, e.g. 'MSS_low_pre'

res    = mdl.Residuals.Raw;       % raw residuals
fitted = mdl.Fitted;              % fitted values
group  = T.isPRT;                 % 0 = control, 1 = PRT

colors = [0 0.447 0.741; 0.85 0.325 0.098];  % MATLAB default blue/orange

%% 1 & 2. Linearity + Homoscedasticity
figure('Color','w','Name','Linearity & Homoscedasticity');

% a) Residuals vs Fitted
subplot(1,2,1); hold on;
gscatter(fitted, res, group, colors, '.', 12);
yline(0,'k--');
xlabel('Fitted values');
ylabel('Raw residuals');
title('Residuals vs Fitted');
grid on; box on;
legend({'Control','PRT'},'Location','best');

% b) Residuals vs MSS
subplot(1,2,2); hold on;
xMSS = T.(mssVarName);
gscatter(xMSS, res, group, colors, '.', 12);
yline(0,'k--');
xlabel(mssVarName,'Interpreter','none');
ylabel('Raw residuals');
title(['Residuals vs ', mssVarName],'Interpreter','none');
grid on; box on;
legend({'Control','PRT'},'Location','best');

%% 3. Normality (QQ-plot)
figure('Color','w','Name','Normality of Residuals');
qqplot(res);
title('QQ-Plot of Residuals');
grid on; box on;

%% 4. Influential points: Cook''s distance & leverage
D   = mdl.Diagnostics.CooksDistance;
lev = mdl.Diagnostics.Leverage;
n   = mdl.NumObservations;
p   = mdl.NumEstimatedCoefficients;

cookCut = 4/n;         % common rule-of-thumb
levCut  = 2*p/n;

figure('Color','w','Name','Influential Points');

% a) Cook's distance vs observation index
subplot(1,2,1); hold on;
stem(1:n, D, 'filled');
yline(cookCut,'r--','4/n','LabelHorizontalAlignment','left');
xlabel('Observation index');
ylabel('Cook''s distance');
title('Cook''s Distance');
grid on; box on;

% b) Leverage vs fitted (or index)
subplot(1,2,2); hold on;
scatter(lev, abs(res), 20, 'k','filled');
xline(levCut,'r--','2p/n');
xlabel('Leverage');
ylabel('|Residual|');
title('Leverage vs |Residual|');
grid on; box on;

% --- Extract Cook's distance ---
D = mdl.Diagnostics.CooksDistance;

% Find the most influential observation
[~, idx_out] = max(D);

% Get the corresponding subject ID
outlier_SubID = predictTbl.SubID(idx_out);

fprintf('Removing influential subject with SubID = %d\n', outlier_SubID);

% --- Remove that subject entirely from the table ---
predictTbl_clean = predictTbl(predictTbl.SubID ~= outlier_SubID, :);
predictTbl = predictTbl_clean;

end
