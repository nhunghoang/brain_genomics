function parcellation_hoa_sh_cer(path_output)

%% get HOA cortical labels

hoa_path = fullfile('/data1/rubinov_lab/bin/fsl/data/atlases/HarvardOxford/');

% get harvard oxford data
hoa_cort_data = double(niftiread(fullfile(hoa_path, 'HarvardOxford-cort-maxprob-thr50-2mm.nii.gz')));
hoa_cort_data(~hoa_cort_data) = nan;
hoa_cort_data = {hoa_cort_data, hoa_cort_data+max(hoa_cort_data(:))};
hoa_cort_data{1}( 1:45 , :, :) = nan;
hoa_cort_data{2}(46:end, :, :) = nan;
hoa_cort_data = max(hoa_cort_data{1}, hoa_cort_data{2});

% get harvard oxford labels
hoa_cort_name_raw = fileread(fullfile(hoa_path, '..', 'HarvardOxford-Cortical.xml'));
name_sta = strfind(hoa_cort_name_raw, '<label');
name_fin = strfind(hoa_cort_name_raw, '</label>');

r = numel(name_sta);
assert(isequal(r, numel(name_fin)))
hoa_cort_name = strings(r, 1);
for i = 1:r
    l = split(hoa_cort_name_raw(name_sta(i):name_fin(i)-1), '>');
    hoa_cort_name(i) = l{2};
end
hoa_cort_name = ["Left "+hoa_cort_name; "Right "+hoa_cort_name];

%% get HOA subcortical labels

% get subcortical data
hoa_subc_data0 = double(niftiread(fullfile(hoa_path, 'HarvardOxford-sub-maxprob-thr50-2mm.nii.gz')));
hoa_subc_data0(~hoa_subc_data0) = nan;

% get subcortical  labels
hoa_subc_name_raw = fileread(fullfile(hoa_path, '..', 'HarvardOxford-Subcortical.xml'));
name_sta = strfind(hoa_subc_name_raw, '<label');
name_fin = strfind(hoa_subc_name_raw, '</label>');

r0 = numel(name_sta);
assert(isequal(r0, numel(name_fin)))
hoa_subc_name0 = strings(r0, 1);
for i = 1:r0
    l = split(hoa_subc_name_raw(name_sta(i):name_fin(i)-1), '>');
    hoa_subc_name0(i) = l{2};
end

% get valid parcels (exclude cortex, white matter, ventricles, brainstem)
idx =       ~contains(hoa_subc_name0, "Cerebral White Matter");
idx = idx & ~contains(hoa_subc_name0, "Lateral Ventricle");
idx = idx & ~contains(hoa_subc_name0, "Cerebral Cortex");
idx = idx & ~contains(hoa_subc_name0, "Brain-Stem");
ldx = find(idx);

% reassign parcel names
r = numel(ldx);
hoa_subc_data = nan(size(hoa_subc_data0));
hoa_subc_name = strings(r, 1);
for i = 1:r
    hoa_subc_data(hoa_subc_data0 == ldx(i)) = i;
    hoa_subc_name(i) = hoa_subc_name0(ldx(i));
end

%% get substantia nigra and hypothalamus

rfl_path = fullfile('/data1/datasets/parcellations/CIT168_Reinf_Learn_v1.1.0');

% load data and resample to 2mm
rfl4_data_1mm_4d = niftiread(fullfile(rfl_path, 'MNI152-FSL/CIT168toMNI152-FSL_prob.nii.gz'));
[x, y, z, l] = size(rfl4_data_1mm_4d);
[Xq, Yq, Zq, Lq] = ndgrid(1.5:2:x, 1.5:2:y, 1.5:2:z, 1:l);
rfl4_data_4d = interpn(rfl4_data_1mm_4d, Xq, Yq, Zq, Lq);

% get parcellations
rfl_data_4d = rfl4_data_4d .* (rfl4_data_4d > 0.5);

rfl_name = readcell(fullfile(rfl_path, 'labels.txt'));
rfl_name = string(rfl_name(:, 2));

% get substantia nigra and hypothalamus

sna_data = any(rfl_data_4d(:, :, :, contains(rfl_name, "SN")), 4);
hth_data = any(rfl_data_4d(:, :, :, contains(rfl_name, "HTH")), 4);

assert(~any(sna_data & hth_data, 'all'))

rfl_data = sna_data + 2*hth_data;
rfl_data(~rfl_data) = nan;

rfl_name = ["SN"; "HTH"];

%% get cerebellar parcellations

cer_path = fullfile('/data1/rubinov_lab/bin/fsl/data/atlases/Cerebellum');

% get parcellations
cer_data_4d = double(niftiread(fullfile(cer_path, 'Cerebellum-MNIfnirt-prob-2mm.nii.gz')));
cer_data_4d = cer_data_4d .* (cer_data_4d > 50);

cer_name_raw = fileread(fullfile(cer_path, '/../Cerebellum_MNIfnirt.xml'));
name_sta = strfind(cer_name_raw, '<label');
name_fin = strfind(cer_name_raw, '</label>');

r = numel(name_sta);
assert(isequal(r, numel(name_fin)))
cer_name = strings(r, 1);
for i = 1:r
    l = split(cer_name_raw(name_sta(i):name_fin(i)-1), '>');
    cer_name(i) = l{2};
end

cer_name_r = contains(cer_name, "Right");
cer_name_l = contains(cer_name, "Left");
cer_name_v = ~(cer_name_l | cer_name_r);

cer_data_r = any(cer_data_4d(:, :, :, cer_name_r), 4);
cer_data_l = any(cer_data_4d(:, :, :, cer_name_l), 4);
cer_data_v = any(cer_data_4d(:, :, :, cer_name_v), 4);

assert(~any(cer_data_r & cer_data_l, 'all'));
assert(~any(cer_data_r & cer_data_v, 'all'));
assert(~any(cer_data_l & cer_data_v, 'all'));

cer_data = cer_data_r + 2 * cer_data_l + 3 * cer_data_v;
cer_data(~cer_data) = nan;

cer_name = ["Right Cerebellum"; "Left Cerebellum"; "Vermis Cerebellum"];

%%

% parc is the full parcellation
parc = 0;
for pi = {hoa_cort_data, hoa_subc_data, rfl_data, cer_data}
    parc = max(parc, pi{1} + max(parc, [], 'all'));
end

% name is the full label vector
name = [hoa_cort_name(:, 1); hoa_subc_name; rfl_name; cer_name];

% run preprocesing
mkdir(path_output)
save(fullfile(path_output, "parcellation"), "parc", "name");
preprocess_hcp_volume(parc, fullfile(path_output, "timeseries"))
