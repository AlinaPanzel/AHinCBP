function [d, lo] = load_data_objects(d, datadir, lo)

%% === Load S1/S2 (patients) and add HC for S1 baseline  ===

% ---------- Patients: Session 1 (ids sub-0*/sub-1*) ----------
clbpS1_t_l_g1 = filenames(fullfile(datadir, 'sub-0*', 'ses-01', 'acute_v2', 'con_0001.nii')); % thumb low
clbpS1_t_l_g2 = filenames(fullfile(datadir, 'sub-1*', 'ses-01', 'acute_v2', 'con_0001.nii'));
g  = vertcat(clbpS1_t_l_g1, clbpS1_t_l_g2);                    % just to extract IDs present at S1
subs  = str2double(extractBetween(g,"sub-","/ses"));

% ---------- Patients: Session 2 ----------
clbpS2_t_l_g1 = filenames(fullfile(datadir, 'sub-0*', 'ses-02', 'acute_v2', 'con_0001.nii')); % thumb low
clbpS2_t_l_g2 = filenames(fullfile(datadir, 'sub-1*', 'ses-02', 'acute_v2', 'con_0001.nii'));
g2 = vertcat(clbpS2_t_l_g1, clbpS2_t_l_g2);
subs2 = str2double(extractBetween(g2,"sub-","/ses"));

% ---------- Healthy Controls: Session 1 ONLY (ids sub-9*) ----------
hcS1_t_l = filenames(fullfile(datadir, 'sub-9*', 'ses-01', 'acute_v2', 'con_0001.nii')); % thumb low
HC_subs  = str2double(extractBetween(hcS1_t_l,"sub-","/ses"));

% Report
fprintf('Found %d patient IDs with S1 and %d with S2; %d HC at S1\n', ...
    numel(subs), numel(subs2), numel(HC_subs));

% All unique IDs (patients only)
all_subs = unique([subs; subs2]);
fprintf('Total unique patient subjects across both sessions: %d\n', numel(all_subs));

% Groups from metadata (patients 1/2/3)
G1 = d.metadata.id(d.metadata.group == 1);
G2 = d.metadata.id(d.metadata.group == 2);
G3 = d.metadata.id(d.metadata.group == 3);

% Intersections: patients per group, per session
[l.ID.G1_S1, ~] = intersect(subs,  G1);
[l.ID.G2_S1, ~] = intersect(subs,  G2);
[l.ID.G3_S1, ~] = intersect(subs,  G3);

[l.ID.G1_S2, ~] = intersect(subs2, G1);
[l.ID.G2_S2, ~] = intersect(subs2, G2);
[l.ID.G3_S2, ~] = intersect(subs2, G3);

% HC IDs (baseline)
l.ID.HC_S1 = HC_subs(:);

fprintf('Group 1: S1 %d, S2 %d\n', numel(l.ID.G1_S1), numel(l.ID.G1_S2));
fprintf('Group 2: S1 %d, S2 %d\n', numel(l.ID.G2_S1), numel(l.ID.G2_S2));
fprintf('Group 3: S1 %d, S2 %d\n', numel(l.ID.G3_S1), numel(l.ID.G3_S2));
fprintf('HC (baseline S1): %d\n', numel(l.ID.HC_S1));

%% ---------- Build file lists ----------

% Session 1 (patients) — THUMB low/high + SOUND low/high
Treatment_S1 = {'G1_S1','G2_S1','G3_S1'};
for j = 1:numel(Treatment_S1)
    for i = 1:numel(l.ID.(Treatment_S1{j}))
        n = l.ID.(Treatment_S1{j})(i);
        str = "sub-" + sprintf('%04d', n);
        l.(Treatment_S1{j}).S1.all_t_l(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0001.nii'));
        l.(Treatment_S1{j}).S1.all_t_h(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0002.nii'));
        l.(Treatment_S1{j}).S1.all_s_l(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0003.nii'));
        l.(Treatment_S1{j}).S1.all_s_h(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0004.nii'));
    end
end

% Session 1 (HC) — THUMB low/high + SOUND low/high
for i = 1:numel(l.ID.HC_S1)
    n = l.ID.HC_S1(i);
    str = "sub-" + sprintf('%04d', n);
    l.HC_S1.S1.all_t_l(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0001.nii'));
    l.HC_S1.S1.all_t_h(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0002.nii'));
    l.HC_S1.S1.all_s_l(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0003.nii'));
    l.HC_S1.S1.all_s_h(i,:) = filenames(fullfile(datadir, str, 'ses-01', 'acute_v2', 'con_0004.nii'));
end

% Session 2 (patients) — THUMB low/high + SOUND low/high
Treatment_S2 = {'G1_S2','G2_S2','G3_S2'};
for j = 1:numel(Treatment_S2)
    for i = 1:numel(l.ID.(Treatment_S2{j}))
        n = l.ID.(Treatment_S2{j})(i);
        str = "sub-" + sprintf('%04d', n);
        l.(Treatment_S2{j}).S2.all_t_l(i,:) = filenames(fullfile(datadir, str, 'ses-02', 'acute_v2', 'con_0001.nii'));
        l.(Treatment_S2{j}).S2.all_t_h(i,:) = filenames(fullfile(datadir, str, 'ses-02', 'acute_v2', 'con_0002.nii'));
        l.(Treatment_S2{j}).S2.all_s_l(i,:) = filenames(fullfile(datadir, str, 'ses-02', 'acute_v2', 'con_0003.nii'));
        l.(Treatment_S2{j}).S2.all_s_h(i,:) = filenames(fullfile(datadir, str, 'ses-02', 'acute_v2', 'con_0004.nii'));
    end
end

%% ---------- Load fmri_data objects ----------

field = {'all_t_l','all_t_h', 'all_s_l', 'all_s_h'};
TreatmentsPatients = {'G1','G2','G3'};

% Session 1 (patients)
fprintf('Loading fMRI data for Session 1 patients...\n');
for n = 1:numel(TreatmentsPatients)
    tag = [TreatmentsPatients{n} '_S1'];
    for i = 1:numel(field)
        if isfield(l, tag) && isfield(l.(tag), 'S1') && isfield(l.(tag).S1, field{i})
            lo.(TreatmentsPatients{n}).S1.(field{i}) = fmri_data(l.(tag).S1.(field{i}));
            fprintf('Loaded %s S1 %s: %d subjects\n', TreatmentsPatients{n}, field{i}, size(lo.(TreatmentsPatients{n}).S1.(field{i}).dat,2));
        else
            fprintf('Warning: no data for %s S1 %s\n', TreatmentsPatients{n}, field{i});
        end
    end
end

% Session 1 (HC)
fprintf('Loading fMRI data for Session 1 HC...\n');
for i = 1:numel(field)
    if isfield(l,'HC_S1') && isfield(l.HC_S1,'S1') && isfield(l.HC_S1.S1, field{i})
        lo.HC.S1.(field{i}) = fmri_data(l.HC_S1.S1.(field{i}));
        fprintf('Loaded HC S1 %s: %d subjects\n', field{i}, size(lo.HC.S1.(field{i}).dat,2));
    else
        fprintf('Warning: no data for HC S1 %s\n', field{i});
    end
end

% Session 2 (patients)
fprintf('Loading fMRI data for Session 2 patients...\n');
for n = 1:numel(TreatmentsPatients)
    tag = [TreatmentsPatients{n} '_S2'];
    for i = 1:numel(field)
        if isfield(l, tag) && isfield(l.(tag), 'S2') && isfield(l.(tag).S2, field{i})
            lo.(TreatmentsPatients{n}).S2.(field{i}) = fmri_data(l.(tag).S2.(field{i}));
            fprintf('Loaded %s S2 %s: %d subjects\n', TreatmentsPatients{n}, field{i}, size(lo.(TreatmentsPatients{n}).S2.(field{i}).dat,2));
        else
            fprintf('Warning: no data for %s S2 %s\n', TreatmentsPatients{n}, field{i});
        end
    end
end

%% ---------- Masks & metadata filtering (based on SOUND presence) ----------
is_nonempty = @(x) ~(isempty(x) || (isstring(x) && x=="") || (iscell(x) && (isempty(x{1}) || strcmp(x{1},""))));

% Session 1 masks (patients)
mask.G1.S1 = arrayfun(is_nonempty, l.G1_S1.S1.all_s_l) & arrayfun(is_nonempty, l.G1_S1.S1.all_s_h);
mask.G2.S1 = arrayfun(is_nonempty, l.G2_S1.S1.all_s_l) & arrayfun(is_nonempty, l.G2_S1.S1.all_s_h);
mask.G3.S1 = arrayfun(is_nonempty, l.G3_S1.S1.all_s_l) & arrayfun(is_nonempty, l.G3_S1.S1.all_s_h);

validIDs.G1.S1 = l.ID.G1_S1(mask.G1.S1);
validIDs.G2.S1 = l.ID.G2_S1(mask.G2.S1);
validIDs.G3.S1 = l.ID.G3_S1(mask.G3.S1);

% Session 1 masks (HC)
mask.HC.S1 = arrayfun(is_nonempty, l.HC_S1.S1.all_s_l) & arrayfun(is_nonempty, l.HC_S1.S1.all_s_h);
validIDs.HC.S1 = l.ID.HC_S1(mask.HC.S1);

% Session 2 masks (patients)
mask.G1.S2 = arrayfun(is_nonempty, l.G1_S2.S2.all_s_l) & arrayfun(is_nonempty, l.G1_S2.S2.all_s_h);
mask.G2.S2 = arrayfun(is_nonempty, l.G2_S2.S2.all_s_l) & arrayfun(is_nonempty, l.G2_S2.S2.all_s_h);
mask.G3.S2 = arrayfun(is_nonempty, l.G3_S2.S2.all_s_l) & arrayfun(is_nonempty, l.G3_S2.S2.all_s_h);

validIDs.G1.S2 = l.ID.G1_S2(mask.G1.S2);
validIDs.G2.S2 = l.ID.G2_S2(mask.G2.S2);
validIDs.G3.S2 = l.ID.G3_S2(mask.G3.S2);

% Completers (patients only)
validIDs.G1.S12 = intersect(validIDs.G1.S1, validIDs.G1.S2);
validIDs.G2.S12 = intersect(validIDs.G2.S1, validIDs.G2.S2);
validIDs.G3.S12 = intersect(validIDs.G3.S1, validIDs.G3.S2);

%% ---------- METADATA SLICES (FIXED HC HANDLING) ----------
meta = d.metadata;

% Robust HC detection:
% * numeric group: HC coded as NaN -> use isnan
% * categorical/char fallback: treat missing/undefined/"NaN" as HC
if isnumeric(meta.group)
    isHC = isnan(meta.group);
elseif iscategorical(meta.group)
    isHC = isundefined(meta.group) | (string(meta.group)=="NaN");
else
    gstr = string(meta.group);
    isHC = (gstr=="NaN") | (gstr=="") | (gstr=="<undefined>");
end

% Patients by group code (1/2/3) at each session, filtered by valid IDs
d.G1.S1 = meta(meta.group==1 & meta.time==1 & ismember(meta.id, validIDs.G1.S1), :);
d.G2.S1 = meta(meta.group==2 & meta.time==1 & ismember(meta.id, validIDs.G2.S1), :);
d.G3.S1 = meta(meta.group==3 & meta.time==1 & ismember(meta.id, validIDs.G3.S1), :);

d.G1.S2 = meta(meta.group==1 & meta.time==2 & ismember(meta.id, validIDs.G1.S2), :);
d.G2.S2 = meta(meta.group==2 & meta.time==2 & ismember(meta.id, validIDs.G2.S2), :);
d.G3.S2 = meta(meta.group==3 & meta.time==2 & ismember(meta.id, validIDs.G3.S2), :);

% Healthies (group coded as NaN/missing), baseline only, filtered by valid HC IDs
d.HC.S1 = meta(isHC & meta.time==1 & ismember(meta.id, validIDs.HC.S1), :);

% Debug prints
fprintf('\nBehavior kept where brain files exist (sound low+high required):\n');
fprintf('G1: S1 %d, S2 %d, S12 %d\n', height(d.G1.S1), height(d.G1.S2), numel(validIDs.G1.S12));
fprintf('G2: S1 %d, S2 %d, S12 %d\n', height(d.G2.S1), height(d.G2.S2), numel(validIDs.G2.S12));
fprintf('G3: S1 %d, S2 %d, S12 %d\n', height(d.G3.S1), height(d.G3.S2), numel(validIDs.G3.S12));
fprintf('HC: S1 %d\n', height(d.HC.S1));
fprintf('Data loading complete (SOUND/THUMB lists built; SOUND masks used).\n');

%% ---------- Completers fMRI (patients only) ----------
lo.G1.S1.completers_s_l = fmri_data(l.G1_S1.S1.all_s_l(ismember(l.ID.G1_S1, validIDs.G1.S12)));
lo.G1.S1.completers_s_h = fmri_data(l.G1_S1.S1.all_s_h(ismember(l.ID.G1_S1, validIDs.G1.S12)));
lo.G1.S2.completers_s_l = fmri_data(l.G1_S2.S2.all_s_l(ismember(l.ID.G1_S2, validIDs.G1.S12)));
lo.G1.S2.completers_s_h = fmri_data(l.G1_S2.S2.all_s_h(ismember(l.ID.G1_S2, validIDs.G1.S12)));

lo.G2.S1.completers_s_l = fmri_data(l.G2_S1.S1.all_s_l(ismember(l.ID.G2_S1, validIDs.G2.S12)));
lo.G2.S1.completers_s_h = fmri_data(l.G2_S1.S1.all_s_h(ismember(l.ID.G2_S1, validIDs.G2.S12)));
lo.G2.S2.completers_s_l = fmri_data(l.G2_S2.S2.all_s_l(ismember(l.ID.G2_S2, validIDs.G2.S12)));
lo.G2.S2.completers_s_h = fmri_data(l.G2_S2.S2.all_s_h(ismember(l.ID.G2_S2, validIDs.G2.S12)));

lo.G3.S1.completers_s_l = fmri_data(l.G3_S1.S1.all_s_l(ismember(l.ID.G3_S1, validIDs.G3.S12)));
lo.G3.S1.completers_s_h = fmri_data(l.G3_S1.S1.all_s_h(ismember(l.ID.G3_S1, validIDs.G3.S12)));
lo.G3.S2.completers_s_l = fmri_data(l.G3_S2.S2.all_s_l(ismember(l.ID.G3_S2, validIDs.G3.S12)));
lo.G3.S2.completers_s_h = fmri_data(l.G3_S2.S2.all_s_h(ismember(l.ID.G3_S2, validIDs.G3.S12)));

end
