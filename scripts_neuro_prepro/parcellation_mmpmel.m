
%%% In this script we get 
%21 sub-cortical regions from: 
%   Melbourne atlas + substantia nigra + hypothalamus + cerebellar
% & 
%360 cortical regions from:
%   Glasser2016_HCP_MMP1.0_5 surface-based atlas

% Output: Parc = cell{1,2}
%         Parc_label = [381,1];

%% ----------- get subcortical parcellations :  21 regions ---------------

%% Melburne parcellation: 16 regions
zal_path = fullfile('/data1/datasets/parcellations/Tian2020MSA_v1.1/3T/Subcortex-Only/');
zal_data = niftiread(fullfile(zal_path, 'Tian_Subcortex_S1_3T.nii.gz'));
zal_data(~zal_data) = nan;

zal_name = string(readcell(fullfile(zal_path, 'Tian_Subcortex_S1_3T_label.txt')));

%% get substantia nigra and hypothalamus : 2 regions

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

%% get cerebellar parcellations : 3 regions

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

%% put subcortical parcellations together

% parc is the full parcellation
Parc = 0;
for pi = {zal_data, rfl_data, cer_data}
    Parc = max(Parc, pi{1} + max(Parc, [], 'all'));
end

%volume data [91*109*91] for subcortical
Parc = double(Parc); 

% subcortical labels
subcortical_name = [zal_name; rfl_name; cer_name];

%% ------------- get cortical parcellations: 360 regions ----------------

% Glasser2016_HCP_MMP1.0_5 surface-based atlas

addpath /data1/rubinov_lab/bin/fieldtrip-20210629
addpath /data1/rubinov_lab/bin/fieldtrip-20210629/fileio
addpath /data1/rubinov_lab/bin/fieldtrip-20210629/utilities

%Glasser_HCP parcellations
cift = ft_read_cifti(fullfile('/data1/rubinov_lab/datasets/human_parcellations/', ...
    'Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'));

% cortical labels
MMP_name = cift.indexmaxlabel;

% Parc : cell{21 Subcortical parcs, 360 Corticall parcs};  Pmax=[21,360];
Parc = {Parc(:), [cift.indexmax(cift.brainstructure==1); cift.indexmax(cift.brainstructure==2)]};
Parc_label = [zal_name; rfl_name; cer_name; MMP_name];

% Ix_parc = cellfun(@(x) find(x > 0), Parc, 'uniformoutput', false);
% N_voxels = cellfun(@nnz, Ix_parc);
% Pmax = cellfun(@(x) max(x(:)), Parc);

%% run preprocesing_hcp_mmpmel

preprocessing_hcp_mmpmel(Parc, Parc_label);



