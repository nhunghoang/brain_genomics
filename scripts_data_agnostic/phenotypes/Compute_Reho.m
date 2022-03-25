

function [Reho_parcelvalues] = Compute_Reho(project, regions, path_output)

%-----------------------------------------------------------------
% example: Reho_parcelvalues = Compute_Reho("hcp", "all", '/data1/rubinov_lab/Neda/ReHo')

%   function summary:
%     this script load preprocessed rs-fmri timeseries as input: Vp_clean
%       >>>>> processing steps
%	     >Extracting ROI timeseries 
%	     >Calculate correlation of timeseries in each ROI: [#voxs, #voxes]
%	     >Get upper/lower triangle values of the symmetric corr matrix
%	     >Calculate mean as the regional homogeneity value

%
%   function inputs:
%     project:      label: human connectome project("hcp") or UK-biobank("ukb")
%     regions:      interested in specific brain regions in brain_genomic project("gene") or 
%                   all brain regions in hoacer_sn_hth atlas("all")
%     path_output:  path to store output for example: path_output = '/data1/rubinov_lab/Neda/Reho'
%
%   function outouts:
%     Reho_parcelvalues :  contain Reho result matrix for defined brain regions :[# regions, # subjects]
%
%
%  The result matrices for "HCP" and "UKB" are saved in the following paths:
%---------UKB------------- 
%  path_output = '/data1/rubinov_lab/Neda/ReHo_UKB/'
%   >> Reho_parcelled.mat; 

%---------HCP-----------
%  path_output = '/data1/rubinov_lab/Neda/ReHo/'
%   >> Reho_parcelled.mat;
%
%*Note: the calculation pipeline for both cohorts is the same!
%
% Nhung Hoang, Neda Sardaripour, 2021 
%-----------------------------------------------------------------
%-----------------------------------------------------------------

atlas = 'hoacer_sn_hth';

% set paths

path_regs = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/naming-115.txt']; 
path_parc = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/parcellation.mat'];
path_mask = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '_voxels/mask.mat'];
display(path_output)

%------------------------------------------------------------------------------------------       
%>>>>> switch between two projects: 1. Human Connectome Project(hcp) 2.UKBiobank(ukb) <<<<<%
%------------------------------------------------------------------------------------------
switch project
    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>Human Connectome Project (HCP)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------
    case "hcp"

        disp(['Reho calculation in HCP dataset has started...']);

        %load rs-fmri timeseries
        path_input = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '_voxels/timeseries'];
        %path_output = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/phenotypes/reho_by_subj'];

        % get list of subjects
        path_scan = ['/data1/rubinov_lab/brain_genomics/data_HCP/' atlas '/timeseries_order.hdf5'];
        % get subject list, parcellation, and brain mask
        subj_list = h5read(path_scan, '/subjects');
        parc = load(path_parc).parc;
        parc_label = load(path_parc).name;
        mask = load(path_mask).mask;
        parc_masked = parc(mask);

        % get voxel indices per brain regions
        reg_voxels = {};
        switch regions

            %----------------------------------------------------------------------
            %   brain region with gene expression information! (18 brain regions)
            %----------------------------------------------------------------------
            case "gene"
                display(['>>> PrediXcan regions <<<'])
                reg_data = tdfread(path_regs);
                for r = 1:numel(reg_data.INDEX)
                    r_idx = reg_data.INDEX(r) + 1;
                    r_nam = strtrim(replace(reg_data.ABBREV(r,:), '-', '')); % remove hashes and whitespaces in name
                    reg_voxels.(r_nam) = find(parc_masked == r_idx);
                end
                %----------------------------------------------------------------------
                %   All brain regions in hoace_sn_hth atlas (115 brain regions)
                %----------------------------------------------------------------------
            case "all"
                display(['>>> 115 brain regions <<<'])
                for r = 1:max(parc(:))
                    r_nam = strtrim(replace(convertStringsToChars(parc_label(r,:)), '-', '')); % remove hashes and whitespaces in name
                    r_nam = r_nam(find(~isspace(r_nam)));
                    r_nam = regexprep(r_nam,'[,]','_'); %replace ',' with '_'
                    r_nam = regexprep(r_nam,')','_');
                    r_nam = regexprep(r_nam,'(','_');
                    r_nam = erase(r_nam,"'");
                    reg_voxels.(r_nam) = find(parc_masked == r);
                end
        end


        display(['loop over subjects'])
        switch regions
            case "gene"
                Reho_allsubjs = nan(numel(reg_data.INDEX),length(subj_list)); % (# regs * #subjs)
            case "all"
                Reho_allsubjs = nan(numel(fieldnames(reg_voxels)),length(subj_list));
        end

        for i = 1:length(subj_list)
            subj = string(subj_list(i));
            disp([num2str(i), subj]);

            filename = fullfile(path_output, subj + ".mat");
            handle = matfile(filename, 'writable', true);
            Vp_clean = load(fullfile(path_input, subj + ".mat")).Vp_clean;

            % loop over subject scans
            scans = h5read(path_scan, ['/' + subj]);
            switch regions
                case "gene"
                    rehos = nan(length(scans), numel(reg_data.INDEX)); % (# scans * # regs)
                case "all"
                    rehos = nan(length(scans), numel(fieldnames(reg_voxels))); 
            end

            for j = 1:length(scans)
                scan = scans(j) + 1;
                data = Vp_clean{scan};

                switch regions

                    case "gene"

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

                    case "all"
                        
                        % loop over brain regions
                        for r = 1:numel(fieldnames(reg_voxels))
                            r_nam = fieldnames(reg_voxels);
                            r_vox = reg_voxels.(r_nam{r});
                            % compute ReHo
                            r_data = data(r_vox,:);
                            rho = corr(r_data.');
                            tril_idx = tril(true(size(rho)), -1);
                            rehos(j,r) = mean(rho(tril_idx), 'omitnan');
                        end
                end
            end

            % compute avg ReHo (across scans, per region)
            mean_reho = mean(rehos, 1, 'omitnan');
            handle.reho = mean_reho;
            disp([subj+' finished'])

            %store all subject data in a matrix [#regions, #subjects]
            Reho_allsubjs(:,i)= mean_reho';


        end %over subjs
        display(['HCP Reho calculation has finished!'])

        %>>>>>>>output variable of function<<<<<<<
        Reho_parcelvalues = Reho_allsubjs;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        %--------------------------------------------------------------------------------------
        %       >>>>>>>>>>>>>>>>>>>>UK Biobank (UKB)<<<<<<<<<<<<<<<<<<<<<
        %--------------------------------------------------------------------------------------
    case "ukb"

        disp(['Reho calculation in ukb dataset has started...']);

        path_input = ['/data1/rubinov_lab/brain_genomics/data_UKB/' atlas '_voxels/timeseries'];
        %path_output = ['/data1/rubinov_lab/brain_genomics/data_UKB/' atlas '/phenotypes/reho_by_subj'];

        % get subject list, parcellation, and brain mask  
        path_scan = ['/data1/rubinov_lab/brain_genomics/data_UKB/' atlas '/timeseries_order.hdf5'];
        subj_list = h5read(path_scan, '/subjects');
        parc = load(path_parc).parc;
        parc_label = load(path_parc).name;
        mask = load(path_mask).mask; 
        parc_masked = parc(mask);

        % get voxel indices per brain regions
        reg_voxels = {};
        switch regions
            %----------------------------------------------------------------------
            %   brain region with gene expression information! (18 brain regions)
            %----------------------------------------------------------------------
            % get voxel indices per PrediXcan region 

            case "gene"
                display(['>>> PrediXcan regions <<<'])
                reg_data = tdfread(path_regs);
                for r = 1:numel(reg_data.INDEX)
                    r_idx = reg_data.INDEX(r) + 1;
                    r_nam = strtrim(replace(reg_data.ABBREV(r,:), '-', '')); % remove hashes and whitespaces in name
                    reg_voxels.(r_nam) = find(parc_masked == r_idx);
                end
                %----------------------------------------------------------------------
                %   All brain regions in hoace_sn_hth atlas (115 brain regions)
                %----------------------------------------------------------------------
            case "all"
                display(['>>> 115 brain regions <<<'])
                for r = 1:max(parc(:))
                    r_nam = strtrim(replace(convertStringsToChars(parc_label(r,:)), '-', '')); % remove hashes and whitespaces in name
                    r_nam = r_nam(find(~isspace(r_nam)));
                    r_nam = regexprep(r_nam,'[,]','_'); %replace ',' with '_'
                    r_nam = regexprep(r_nam,')','_');
                    r_nam = regexprep(r_nam,'(','_');
                    r_nam = erase(r_nam,"'");
                    reg_voxels.(r_nam) = find(parc_masked == r);
                end
        end

        display(['loop over subjects'])
        switch regions
            case "gene"
                Reho_allsubjs = nan(numel(reg_data.INDEX),length(subj_list)); % (# regs * #subjs)
            case "all"
                Reho_allsubjs = nan(numel(fieldnames(reg_voxels)),length(subj_list));
        end
        for i = 1:length(subj_list)
            subj = string(subj_list(i));
            filename = fullfile(path_output, subj + ".mat");
            handle = matfile(filename, 'writable', true);
            Vp_clean = load(fullfile(path_input, subj + ".mat")).Vp_clean;

            % read subject scan  --> one scan in ukb
            data = Vp_clean{1};

                switch regions

                    case "gene"

                        % loop over PrediXcan regions
                        rehos = nan(numel(reg_data.INDEX), 1); % size: (# regs, 1)  

                        for r = 1:numel(reg_data.INDEX)
                            r_nam = strtrim(replace(reg_data.ABBREV(r,:), '-', ''));
                            r_vox = reg_voxels.(r_nam);
                            % compute ReHo
                            r_data = data(r_vox,:);
                            rho = corr(r_data.');
                            tril_idx = tril(true(size(rho)), -1);
                            rehos(r) = mean(rho(tril_idx), 'omitnan');
                        end


                    case "all"

                        % loop over brain regions
                        rehos = nan(numel(fieldnames(reg_voxels)), 1); 

                        for r = 1:numel(fieldnames(reg_voxels))
                            r_nam = fieldnames(reg_voxels);
                            r_vox = reg_voxels.(r_nam{r});
                            % compute ReHo
                            r_data = data(r_vox,:);
                            rho = corr(r_data.');
                            tril_idx = tril(true(size(rho)), -1);
                            rehos(r) = mean(rho(tril_idx), 'omitnan');
                        end
                end

            % write regional ReHo values 
            handle.reho = rehos; 
            disp([subj+' finished'])

            %store all subject data in a matrix [#regions, #subjects]
            Reho_allsubjs(:,i)= rehos';
        end %over subjs
        display(['UKB Reho calculation has finished!'])

        %>>>>>>>output variable of function<<<<<<<
        Reho_parcelvalues = Reho_allsubjs;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
end

        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/Reho

        save(fullfile(path_output,'Reho_parcelled.mat'),'Reho_allsubjs','-v7.3') % Reho matrix 
        disp(['parcel based result has saved in output_path']);


