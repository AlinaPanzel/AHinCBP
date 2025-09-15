%% ===== SOUND =====
disp('==== SOUND TTEST ====');
sound_baseline_roi = neural_baseline_roi_clean(neural_baseline_roi_clean.modality=="Sound",:);
T = sound_baseline_roi;
T.GroupBin = categorical(ismember(T.group,[1 2 3]), [0 1], {'HC','CBP'});

print_roi_block(T);

%% ===== PRESSURE =====
disp('==== PRESSURE TTEST ====');
pressure_baseline_roi = neural_baseline_roi_clean(neural_baseline_roi_clean.modality=="Pressure",:);
T = pressure_baseline_roi;
T.GroupBin = categorical(ismember(T.group,[1 2 3]), [0 1], {'HC','CBP'});   % (was missing)

print_roi_block(T);

%% ---------- helper ----------
function print_roi_block(T)
ROIs = categories(categorical(T.measure));
ints = categories(categorical(T.intensity));

for ii = 1:numel(ints)
    thisInt = ints{ii};
    fprintf('\n--- Intensity: %s ---\n', string(thisInt));

    % header
    fprintf('%-16s %6s %6s %10s %10s %8s %10s %8s %4s\n', ...
        'ROI','nHC','nCBP','meanHC','meanCBP','t','p','q(BH)','sig');
    fprintf('%s\n', repmat('-',1,88));

    p = nan(numel(ROIs),1);
    rows = strings(numel(ROIs),1);

    for r = 1:numel(ROIs)
        Tk = T(T.measure==ROIs{r} & T.intensity==thisInt,:);

        x = (Tk.value(Tk.GroupBin=="HC"));
        y = (Tk.value(Tk.GroupBin=="CBP"));

        if isempty(x) || isempty(y)
            p(r) = NaN;
            rows(r) = sprintf('%-16s %6d %6d %10s %10s %8s %10s %8s %4s', ...
                ROIs{r}, numel(x), numel(y), 'NA','NA','NA','NA','NA','');
            continue
        end

        [~,p(r),~,st] = ttest2(x,y);
        rows(r) = sprintf('%-16s %6d %6d %10.2f %10.2f %8.2f %10.4g %8s %4s', ...
            ROIs{r}, numel(x), numel(y), mean(x,'omitnan'), mean(y,'omitnan'), ...
            st.tstat, p(r), '----', '');  % q filled below
    end

    % FDR (Benjamini–Hochberg)
    [~,~,~,q] = fdr_bh(p);

    % print with q + stars
    for r = 1:numel(ROIs)
        qtxt = ifelse(isnan(q(r)),'NA',sprintf('%.4g',q(r)));
        sig  = stars(q(r));
        % replace the placeholder '----' with q
        line = replace(rows(r),"    ----", sprintf('%8s', qtxt));
        fprintf('%s %4s\n', line, sig);
    end
end
end

function s = stars(q)
if isnan(q); s = ''; return; end
if q < 0.001, s = '***';
elseif q < 0.01, s = '**';
elseif q < 0.05, s = '*';
else, s = '';
end
end

function out = ifelse(tf, a, b)
if tf, out = a; else, out = b; end
end
