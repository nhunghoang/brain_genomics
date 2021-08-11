function compute_reho(atlas)

% example inputs 
% atlas = 'hoacer_sn_hth';

% set paths
path_scan = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/timeseries_order.hdf5'];
path_regs = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/naming-115.txt']; 
path_parc = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/parcellation.mat'];
path_mask = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '_voxels/mask.mat'];

path_input = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '_voxels/timeseries'];
path_output = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/phenotypes/reho_by_subj'];
%mkdir(path_output);

% get subject list, parcellation, and brain mask  
subj_list = h5read(path_scan, '/subjects'); 
parc = load(path_parc).parc;
mask = load(path_mask).mask; 

% get voxel indices per PrediXcan region 
reg_data = tdfread(path_regs);
parc_masked = parc(mask); 
reg_voxels = {}; 
for r = 1:numel(reg_data.INDEX)
    r_idx = reg_data.INDEX(r) + 1;
    r_nam = strtrim(replace(reg_data.ABBREV(r,:), '-', '')); % remove hashes and whitespaces in name
    reg_voxels.(r_nam) = find(parc_masked == r_idx); 
end

% loop over subjects 
parfor i = 1:length(subj_list)
    subj = string(subj_list(i)); 
    filename = fullfile(path_output, subj + ".mat"); 
    handle = matfile(filename, 'writable', true);  
    Vp_clean = load(fullfile(path_input, subj + ".mat")).Vp_clean;

    % loop over subject scans 
    scans = h5read(path_scan, ['/' + subj]);
    rehos = nan(length(scans), numel(reg_data.INDEX)); % (# scans * # regs)  
    for j = 1:length(scans)
        scan = scans(j) + 1; 
        data = Vp_clean{scan};
        
        % loop over PrediXcan regions 
        for r = 1:numel(reg_data.INDEX)
            r_nam = strtrim(replace(reg_data.ABBREV(r,:), '-', '')); 
            r_vox = reg_voxels.(r_nam); 

            % compute ReHo 
            r_data = data(r_vox,:); 
            rho = corr(r_data.');
            tril_idx = tril(true(size(rho)), -1);
            rehos(j,r) = mean(rho(tril_idx), 'omitnan'); 
        end
    end

    % compute avg ReHo (across scans, per region) 
    mean_reho = mean(rehos, 1, 'omitnan');
    handle.reho = mean_reho; 
    disp(i)
end 
