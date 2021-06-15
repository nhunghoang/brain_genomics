% paths  
path_input = '/data1/datasets/hcp';
path_output = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/regional_homogeneity2.mat';
subj_output = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/ReHo_by_subj2.mat';

% load parcellation
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/parc-121.mat');
parc = parc_file.parc;
parc = parc(:);

% load region indices of interest 
data = tdfread('/data1/rubinov_lab/brain_genomics/data_HCP/naming-121.txt');
regs = strings(18,1);
for r = 1:18
    reg_str = convertCharsToStrings(data.ABBREV(r,:));
    regs(r) = strrep(strrep(reg_str, '-', '_'), ' ', '');
end
idxs = int8(data.INDEX) + 1;

% get file info of valid subject scans 
ts_order = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5';
subjs = h5read(ts_order, '/subjects');
subj_scans = containers.Map('KeyType', 'char', 'ValueType', 'any');
for s = 1:length(subjs)
    subj = int2str(subjs(s));
    subj_scans(subj) = h5read(ts_order, ['/' subj]);
end 

% get scans of interest
scans = cell2table({ ...
    'rfMRI_REST1_LR', 0.72;
    'rfMRI_REST1_RL', 0.72;
    'rfMRI_REST2_LR', 0.72;
    'rfMRI_REST2_RL', 0.72}, ...
    'variablenames', {'Name', 'TR'});

% loop through subjects 
delete(subj_output);
sfile = matfile(subj_output);
subj_reg_rehos = zeros(length(subjs), length(regs));
for s = 1:length(subjs)
    subj = int2str(subjs(s));
    disp([num2str(s), ': ', subj]);

    % loop through valid scans 
    s_scans = subj_scans(subj) + 1;
    subj_rehos = zeros(length(s_scans), length(regs));
    for ss = 1:length(s_scans)
        scan_type = scans.Name{s_scans(ss)};

        % take the timeseries 
        scan_path = fullfile(path_input, subj, 'MNINonLinear', 'Results', scan_type);
        scan_name = [scan_type '_hp2000_clean.nii.gz'];
        V = niftiread(fullfile(scan_path, scan_name));
        t = size(V, 4); 
        V = double(reshape(V, [], t));

        % loop through regions
        for r = 1:length(regs)
            r_idx = idxs(r);
            reg_sigs = V(parc==r_idx, :);
            
            % ReHo := avg correlation between voxels 
            % nb: NaN is caused by voxels with constant signals 
            vox_corrs = corr(reg_sigs.'); 
            this_idx = tril(true(size(vox_corrs)), -1);
            this_reho = mean(vox_corrs(this_idx), 'omitnan');
            subj_rehos(ss,r) = this_reho;
        end
    end
    
    % average the subject's ReHos across their scans
    avg_subj_rehos = mean(subj_rehos, 1);
    subj_reg_rehos(s,:) = avg_subj_rehos;
    disp(avg_subj_rehos);
    
    % cautionary save
    sfile.(['s_' subj]) = avg_subj_rehos;
end

% loop through regions, save subject array of rehos 
delete(path_output);
rfile = matfile(path_output);
for r = 1:length(regs)
    reg = regs(r);
    rfile.(reg) = subj_reg_rehos(:,r);
end

