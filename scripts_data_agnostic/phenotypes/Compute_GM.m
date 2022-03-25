

function [GM_voxelvalues, GM_parcelvalues] = Compute_GM(project, path_output)

%-----------------------------------------------------------------
%
%   function summary:
%     this script load T1-wieghted MRI 3D-data as input
%       >>>>> processing steps:
%        >Extracting brain regions using BET algorithm 
%	     >Segmenting WM,GM, and CSF using FAST algorithm
%	     >Resampling data to MNI152_2mm space (91*109*91)
%	     >Load data containing GM volume
%
%   function inputs:
%     project:      label: human connectome project("hcp") or UK-biobank("ukb")
%     path_output:  path to store output for example: path_output = '/data1/rubinov_lab/Neda/GM'
%
%   function outouts:
%     GM_parcelvalues :  contain Gray matter result matrix for 115 brain regions based on "hoacer_sn_hth" parcellation file :[115, number of subjects]
%     GM_voxelvalues  :  contain Gray matter result matrix for all brain voxels: [number of voxels, number of subjects]
%
%
%  The result matrices for "HCP" and "UKB" are saved in the following paths:
%---------UKB------------- 
%  path_output = '/data1/rubinov_lab/Neda/GM_UKB/'
%   >> GM_all_voxels.mat; 
%   >> GM_vol_hoacer_sn_hth.mat;

%---------HCP-----------
%  path_output = '/data1/rubinov_lab/Neda/GM/'
%   >> GM_all_voxels.mat; 
%   >> GM_vol_hoacer_sn_hth.mat;
%
%*Note: the calculation pipeline for both cohorts is the same!
%*Note: commented lines have already ran and the output files have generated using FSL. Not necessary to run them egain
%
% Neda Sardaripour, 2021 
%-----------------------------------------------------------------
%-----------------------------------------------------------------

%----- load parcellation: hoacer_sn_hth (115 brain regions) ------
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;
        

%------------------------------------------------------------------------------------------       
%>>>>> switch between two projects: 1. Human Connectome Project(hcp) 2.UKBiobank(ukb) <<<<<%
%------------------------------------------------------------------------------------------
switch project
    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>Human Connectome Project (HCP)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------
    case "hcp"

        disp(['Gray matter calculation in HCP dataset has started...']);

        %load T1w volume (bias corrected)
        path_input = '/data1/datasets/hcp'

        % get list of subjects
        ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5';%890
        subjs = h5read(ts_order, '/subjects');


        % loop through hcp subjects
        GM_vol_allvoxels = nan(length(parc(:)), length(subjs));
        N=0;
        for s = 1:length(subjs)

            subj = num2str(subjs(s));
            disp([num2str(s), ': ', subj]);

            scan_path = fullfile(path_input, subj, 'T1w');
            scan_name = 'T1w_acpc_dc_restore_brain.nii.gz';
            %T1w_brain = niftiread(fullfile(scan_path, scan_name));
            %%----------------------------------------
            %           brain extraction (BET)
            %%------------------------------------------
            %bet = "/data1/rubinov_lab/bin/fsl/bin/bet";
            %command = char([bet+" "+fullfile(scan_path, 'T1w_acpc_dc_restore_brain.nii.gz')+...
            %              " "+ fullfile(T1w_path,'T1w_acpc_dc_restore_brain.nii.gz')+...
            %              " -m"+...
            %              " -f .3"]);
            %          [status,cmdout] = system(command);
            %%------------------------------------------
            %       Segmentation using FAST in FSL
            %%------------------------------------------
            %fast = "/data1/rubinov_lab/bin/fsl/bin/fast";
             write_path = fullfile('/data1/rubinov_lab/datasets/hcp/data', subj);
            %
            %command = char([fast+" -o "+fullfile(write_path,'T1w_acpc_dc_restore_brain')+" "+...
            %               fullfile(scan_path, 'T1w_acpc_dc_restore_brain.nii.gz')]);
            %          [status,cmdout] = system(command);

            %--------------------------------------------------------------
            %       >> output: pve_0:CSF, pve_1:GM, pve_2:WM
            %--------%loading pve_1 gives us Gray matter value-->[0,1]

            if size(dir(fullfile(path_input, subj,'T1w')),1)>3
                %---------------------T1w_seg_GM   :   260x311x260 -------------- ----
                %T1w_seg_GM = niftiread(fullfile(write_path,'T1w_acpc_dc_restore_brain_pve_1.nii.gz'));
                %---------------------------------------------------------
                %           resampling data to MNI152_2mm space
                %%--------------------------------------------------------
                %MNI_standard_path = '/data1/rubinov_lab/bin/fsl/data/standard';

                %------ first compute transform matrix using flirt in FSL
                %flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
                %
                %command = char([flirt+" -in "+fullfile(scan_path, scan_name)+...
                %              " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                %              " -omat "+ fullfile(write_path,'transform_T1w_to_MNI2mm.mat')]);
                %          [status,cmdout] = system(command);
                %
                %-------% apply transform
                %command = char([flirt+" -in "+fullfile(scan_path, 'T1_brain_pve_1.nii.gz')+...
                %              " -ref "+ fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                %              " -omat "+ fullfile(write_path, 'transform_T1_to_MNI2mm.mat')+...
                %              " -out "+ fullfile(write_path, 'T1_brain_pve1_2mm.nii.gz')]);
                %         [status,cmdout] = system(command);


                %-----------reading data--------------
                %check if data exist for rading
                if exist(fullfile(write_path,'T1w_acpc_dc_restore_brain_pve1_2mm.nii.gz'))

                    GM_vol = double(niftiread(fullfile(write_path,'T1w_acpc_dc_restore_brain_pve1_2mm.nii.gz')));
                    GM_vol = GM_vol(:);
                    GM_vol_allvoxels(:,s) = GM_vol;
    
                else
                    N=N+1;
                    no_data(N,1) = num2cell(s);
                    no_data(N,2) = cellstr(subj);
                    disp([num2str(s), ': ', subj,' no data']);
                end
            end
        end
         display(['HCP GM calculation has finished!'])

        %>>>>>>>output variable of function<<<<<<<
        GM_voxelvalues = GM_vol_allvoxels;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/GM

        %voxel-based GM  values
        save(fullfile(path_output,'GM_allvoxels.mat'),'GM_vol_allvoxels','no_data','-v7.3') % ALFF matrix 
        disp(['voxel based result has saved in output_path']);


        %--------------------------------------------------------------------------------------
        %       >>>>>>>>>>>>>>>>>>>>UK Biobank (UKB)<<<<<<<<<<<<<<<<<<<<<
        %--------------------------------------------------------------------------------------
    case "ukb"

        disp(['ALFF/fALFF calculation in ukb dataset has started...']);

        %load T1w volume
        path_input = fullfile('/data1/datasets/ukbiobank/data/T1')
        GM_path = fullfile('/T1/T1_fast');


        % get list of subjects
        UKB_subjs_order = fullfile('/data1/rubinov_lab/Neda/UKB_subjs_list/');
        subjs = load(fullfile(UKB_subjs_order, 'final_subjs.txt'));

        % loop through hcp subjects
        GM_vol_allvoxels = nan(length(parc(:)), length(subjs));
        N=0;
        for s = 1:length(subjs)

            subj = num2str(subjs(s));
            disp([num2str(s), ': ', subj]);

            %------------segmentation result has provided by ukb-----------
            %------------ > output: pve_0:CSF, pve_1:GM, pve_2:WM
            %------------ > loading pve_1 gives us Gray matter value-->[0,1]

            scan_path = fullfile(path_input, subj, GM_path);

            if size(dir(fullfile(scan_path)),1)>3
                %---------------------T1w_seg_GM   :   164x222x185 -------------- ----
                %T1w_seg_GM = niftiread(fullfile(scan_path,'T1_brain_pve_1.nii.gz'));
                %
                %---------------------------------------------------------
                %           resampling data to MNI152_2mm space
                %%--------------------------------------------------------
                %MNI_standard_path = '/data1/rubinov_lab/bin/fsl/data/standard';
                write_path = fullfile('/data1/rubinov_lab/Neda/GM_UKB',subj);
                %-------first compute transform matrix using flirt in FSL
                %flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
                %
                %command = char([flirt+" -in "+fullfile(scan_path, 'T1_brain_pve_1.nii.gz')+...
                %              " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                %              " -omat "+ fullfile(write_path,'transform_T1w_to_MNI2mm.mat')]);
                %         [status,cmdout] = system(command);
                %
                %------% apply transform
                %command = char([flirt+" -in "+fullfile(scan_path, 'T1_brain_pve_1.nii.gz')+...
                %              " -ref "+ fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                %              " -omat "+ fullfile(write_path, 'transform_T1_to_MNI2mm.mat')+...
                %              " -out "+ fullfile(write_path, 'T1_brain_pve1_2mm.nii.gz')]);
                %          [status,cmdout] = system(command);

                %-----------reading data--------------
                %check if data exist for rading
                if exist(fullfile(write_path,'T1_brain_pve1_2mm.nii.gz'))

                    GM_vol = double(niftiread(fullfile(write_path,'T1_brain_pve1_2mm.nii.gz')));
                    GM_vol = GM_vol(:);
                    GM_vol_allvoxels(:,s) = GM_vol;
                else
                    N=N+1;
                    no_data(N,1) = num2cell(s);
                    no_data(N,2) = cellstr(subj);

                    disp([num2str(s), ': ', subj,' no data']);

                end
            end
        end

        display(['ukb GM calculation has finished!'])

         %>>>>>>>output variable of function<<<<<<<
        GM_voxelvalues = GM_vol_allvoxels;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/GM_UKB

        %voxel-based GM  values
        save(fullfile(path_output,'GM_allvoxels.mat'),'GM_vol_allvoxels','no_data','-v7.3') % ALFF matrix 
        disp(['voxel based result has saved in output_path']);

end

%  >> After calculating and saving voxel based GM in HCP or UKB cohorts:
%  >> calculate values for each brain parcels (hoacer_sn_hth parcellation file)
%
% %% loading results
% %voxel data
% GM_vol_allvoxels = load(fullfile('/data1/rubinov_lab/Neda/GM_UKB/GM_vol_allvoxels).GM_vol_allvoxels;

%--------------------------------------------------------------------------------
% -------parcellation: caluculate parcel based GM values------------------
%--------------------------------------------------------------------------------
disp(['parcellation has started...']);

% get mask indices
idx_parc = find(parc);
nmax = nnz(parc);
mask = parc(idx_parc);

 % parcellation
  for s=1:length(subjs)  
      
        for j = 1:max(parc(:))
            GM_115parc(j,s)= sum(GM_vol_allvoxels(parc(:)==j,s));
        end
  end
  save('GM_vol_hoacer_sn_hth.mat','GM_115parc','no_data','-v7.3')  


        
 %>>>>>>>>output variable of function<<<<<<<<<<<<<<
  GM_parcelvalues = GM_115parc;

 %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%----------save data in path_output----------
%parcel based GM
save(fullfile(path_output, 'GM_hoacer_sn_hth_vox_parcelled.mat'),'GM_parcelvalues','no_data','-v7.3') 
disp(['parcellated result has saved in output path']);

end




