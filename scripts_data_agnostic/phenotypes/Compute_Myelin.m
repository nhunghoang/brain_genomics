

function [Myelin_voxelvalues, Myelin_parcelvalues] = Compute_Myelin(project, path_output)

%note: if want to change the output path, first uncomment the FSL related commands
%-----------------------------------------------------------------
% example: [Myelin_voxelvalues, Myelin_parcelvalues] = Compute_Myelin("hcp", '/data1/rubinov_lab/Neda/Myelin')

%   function summary:
%     this script load T1w and T2w MRI data
%       >>>>> processing steps:
%	     >Load T1 and T2 data 
%	     >Transfer data to MNI-2mm space: [91*019*91]
%	     >Clamping range of numbers: [0,100] to avoid generating Inf
%	     >Calculate myelination map by deviding T1 / T2 vols
%
%   function inputs:
%     project:      label: human connectome project("hcp") or UK-biobank("ukb")
%     path_output:  path to store output for example: path_output = '/data1/rubinov_lab/Neda/Myelin'
%
%   function outouts:
%     Myelin_parcelvalues :  contain Myelination result matrix for 115 brain regions based on "hoacer_sn_hth" parcellation file :[115, number of subjects]
%     Myelin_voxelvalues  :  contain Myelination result matrix for all brain voxels: [number of voxels, number of subjects]
%
%  The result matrices for "HCP" and "UKB" are saved in the following paths:
%---------UKB------------- 
%  path_output = '/data1/rubinov_lab/Neda/Myelin_UKB/'
%   >> Myelin_allvoxles.mat; 
%   >> Myelin_hoacer_sn_hth.mat; 

%---------HCP-----------
%  path_output = '/data1/rubinov_lab/Neda/Myelin/'
%   >> Myelin_allvoxels.mat;
%   >> Myelin_hoacer_sn_hth.mat; 
%
%*Note: commented lines have already ran and the output files have
% generated and exist in above^ paths: (Not needed to run again)
%
%
% Neda Sardaripour, 2021 
%-----------------------------------------------------------------
%-----------------------------------------------------------------


% load parcellation
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;

%MNI_standard_path = fullfile('/data1/rubinov_lab/bin/fsl/data/standard');
display(path_output)


%------------------------------------------------------------------------------------------
%>>>>> switch between two projects: 1. Human Connectome Project(hcp) 2.UKBiobank(ukb) <<<<<%
%------------------------------------------------------------------------------------------
switch project
    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>Human Connectome Project (HCP)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------
    case "hcp"

        disp(['Myelin calculation in HCP dataset has started...']);

        %load T1 & T2
        path_input = '/data1/datasets/hcp';
        %path_output = '/data1/rubinov_lab/Neda/Myelin';

        % get list of subjects
        ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5';%890
        subjs = h5read(ts_order, '/subjects');

        Myelin_allvoxels = nan(length(parc(:)), length(subjs));
        Myelin_hoacer_sn_hth = nan(max(parc), length(subjs));
        N=0;

        %loop through subjects
        for s = 1:length(subjs)
        
            subj = num2str(subjs(s));
            disp([num2str(s), ': ', subj]);
        
            % read data 
                scan_path = fullfile(path_input, subj, 'T1w');
                T1_scan = 'T1w_acpc_dc_restore_brain';
                T2_scan = 'T2w_acpc_dc_restore_brain';
         
                write_path = fullfile(path_output, subj);
            
                if ~exist(write_path)
                    mkdir(write_path)
                end
        
            %if size(dir(fullfile(path_input, T1_path, subj, T1_path)),1)>3 && size(dir(fullfile(path_input, T2_path, subj, T2_path)),1)>3
               
%                 %-------load T1 and T2 data : [260*311*260] -------------
%                 T1_vol = double(niftiread(fullfile(scan_path, T1_scan)));
%                 T2_vol = double(niftiread(fullfile(scan_path, T2_scan)));
%                 
%                 %---------------------------------------------------------
%                 %     resampling data to MNI152_2mm space : [91*109*91]
%                 %%--------------------------------------------------------
%                  flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
%                  %T1w
%                  command = char([flirt+" -in "+fullfile(scan_path, 'T1w_acpc_dc_restore_brain')+...
%                      " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                      " -omat "+ fullfile(write_path, 'transform_T1w_to_MNI2mm.mat')]);
%                  [status,cmdout] = system(command);
%                  
%                  %T2w
%                  command = char([flirt+" -in "+fullfile(scan_path, 'T2w_acpc_dc_restore_brain')+...
%                      " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                      " -omat "+ fullfile(write_path, 'transform_T2w_to_MNI2mm.mat')]);
%                  [status,cmdout] = system(command);
%                    
%                    
%                 %apply the same transform matrix to input
%                 %T1w
%                 command = char([flirt+" -in "+fullfile(scan_path, 'T1w_acpc_dc_restore_brain')+...
%                 " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                 " -out "+fullfile(write_path, 'T1w_acpc_dc_restore_brain_MNI2mm')+...
%                 " -init "+fullfile(write_path, 'transform_T1w_to_MNI2mm.mat')+...
%                 " -applyxfm"]);
%                 [status,cmdout] = system(command); disp(cmdout)
%                 
%                 %T2w
%                 command = char([flirt+" -in "+fullfile(scan_path, 'T2w_acpc_dc_restore_brain')+...
%                 " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                 " -out "+fullfile(write_path, 'T2w_acpc_dc_restore_brain_MNI2mm')+...
%                 " -init "+fullfile(write_path, 'transform_T2w_to_MNI2mm.mat')+...
%                 " -applyxfm"]);
%                 [status,cmdout] = system(command); disp(cmdout)
                        
%               %-------load T1 and T2 data : [91*109*91] -------------
%                 ***Note: T1w_acpc_dc_restore_brain_MNI2mm already exist in
%                   path_output = 'data1/rubinov_lab/Neda/Myelin'

                T1_vol = double(niftiread(fullfile(write_path, 'T1w_acpc_dc_restore_brain_MNI2mm')));               
                T2_vol = double(niftiread(fullfile(write_path, 'T2w_acpc_dc_restore_brain_MNI2mm')));
                
                %for k=1:91; imshowpair(map(:,:,k),MNI(:,:,k));pause;end
        
                %-------clamping numbers: [0,100]---------
                % set all negative numbers to 0
                
                T1_vol = (T1_vol)/prctile(nonzeros(T1_vol), 99);
                T2_vol = (T2_vol)/prctile(nonzeros(T2_vol), 99);
                idx = (T2_vol > 0.1);
                
               %----------------------------------------
                %              Myelin map
               %----------------------------------------
                map = zeros(size(T1_vol));
                map(idx) = double(T1_vol(idx)./T2_vol(idx));
                
                Myelin_allvoxels(:,s) = map(:);         
            
        end %through subjects

        display(['HCP Myelin calculation has finished!'])

        %>>>>>>>output variable of function<<<<<<<
        Myelin_voxelvalues = Myelin_allvoxels;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/Myelin

        %voxel-based Myelin  values
        save(fullfile(path_output,'Myelin_allvoxels.mat'),'Myelin_allvoxels','-v7.3') % Myelin matrix 
        disp(['voxel based result has saved in output_path']);

        %--------------------------------------------------------------------------------------
        %       >>>>>>>>>>>>>>>>>>>>UK Biobank (UKB)<<<<<<<<<<<<<<<<<<<<<
        %--------------------------------------------------------------------------------------
    case "ukb"

        disp(['Myelin calculation in ukb dataset has started...']);

        
        path_input = '/data1/datasets/ukbiobank/data';
        T1_path = '/T1';
        T2_path = '/T2_FLAIR';
        %path_output = '/data1/rubinov_lab/Neda/Myelin_UKB';


        % get file info of valid subject scans 
        UKB_subjs_order = fullfile('/data1/rubinov_lab/Neda/UKB_subjs_list/');
        subjs = load(fullfile(UKB_subjs_order, 'final_subjs.txt'));
                
        Myelin_allvoxels = nan(length(parc(:)), length(subjs));
        Myelin_hoacer_sn_hth = nan(max(parc), length(subjs));
        N=0;

        disp(['loop through subjects']);
        for s = 1:length(subjs)
        
            subj = num2str(subjs(s));
            disp([num2str(s), ': ', subj]);
            
            write_path = fullfile(path_output,subj);
            
            if ~exist(write_path)
                mkdir(write_path)
            end

             %check if data exist
            if size(dir(fullfile(path_input, T1_path, subj, T1_path)),1)>3 && size(dir(fullfile(path_input, T2_path, subj, T2_path)),1)>3

                %                % load T1 and T2 data 
%                 T1_scan_path = fullfile(path_input, T1_path, subj, T1_path);
%                 T2_scan_path = fullfile(path_input, T2_path, subj, T2_path);
%                 
%                 %-----------resampling data to 2mm MNI space using FSL command---------------- 
%                 % first compute transform matrix using flirt in FSL
%                  flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
%                  %T1w
%                  command = char([flirt+" -in "+fullfile(T1_scan_path, 'T1_unbiased_brain.nii.gz')+...
%                      " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                      " -omat "+ fullfile(write_path, 'transform_T1w_to_MNI2mm.mat')]);
%                  [status,cmdout] = system(command);
%                  
%                  %T2w
%                  command = char([flirt+" -in "+fullfile(T2_scan_path, 'T2_FLAIR_unbiased_brain.nii.gz')+...
%                      " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                      " -omat "+ fullfile(write_path, 'transform_T2w_to_MNI2mm.mat')]);
%                  [status,cmdout] = system(command);
%                    
%                    
%                 %apply the same transform matrix to input
%                 %T1w
%                 command = char([flirt+" -in "+fullfile(T1_scan_path, 'T1_unbiased_brain.nii.gz')+...
%                 " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                 " -out "+fullfile(write_path, 'T1_unbiased_brain_MNI2mm.nii')+...
%                 " -init "+fullfile(write_path, 'transform_T1w_to_MNI2mm.mat')+...
%                 " -applyxfm"]);
%                 [status,cmdout] = system(command); disp(cmdout)
%                 
%                 %T2w
%                 command = char([flirt+" -in "+fullfile(T2_scan_path, 'T2_FLAIR_unbiased_brain.nii.gz')+...
%                 " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%                 " -out "+fullfile(write_path, 'T2_FLAIR_unbiased_brain_MNI2mm.nii')+...
%                 " -init "+fullfile(write_path, 'transform_T2w_to_MNI2mm.mat')+...
%                 " -applyxfm"]);
%                 [status,cmdout] = system(command); disp(cmdout)
%                         
                scan_name = 'T1_unbiased_brain_MNI2mm.nii.gz';
                T1_vol = double(niftiread(fullfile(write_path, scan_name)));
                
                scan_name_flair = 'T2_FLAIR_unbiased_brain_MNI2mm.nii.gz';
                T2_vol = double(niftiread(fullfile(write_path, scan_name_flair)));
                % for k=1:150; imshowpair(T1_vol(:,:,k),T2_vol(:,:,k));pause;end

                %------clamping numbers: [0,100]
                % set all negative numbers to 0
                
                T1_vol = (T1_vol)/prctile(nonzeros(T1_vol), 99);
                T2_vol = (T2_vol)/prctile(nonzeros(T2_vol), 99);
                
                idx = (T2_vol > 0.1);
                
               %----------------------------------------
                %       Myelin map
               %----------------------------------------
                map = zeros(size(T1_vol));
                map(idx) = double(T1_vol(idx)./T2_vol(idx));
                       
            else 
                    N=N+1;
                    no_data(N,1) = num2cell(s);
                    no_data(N,2) = cellstr(subj);
                    disp([num2str(s), ': ', subj,' no data']); %N=176
            end

           Myelin_allvoxels(:,s) = map(:);  

        end %through subjects

        display(['UKB Myelin calculation has finished!'])

        %>>>>>>>output variable of function<<<<<<<
        Myelin_voxelvalues = Myelin_allvoxels;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/Myelin_UKB

        %voxel-based Myelin  values
        save(fullfile(path_output,'Myelin_allvoxels.mat'),'Myelin_allvoxels','no_data','-v7.3'); % Myelin matrix
        disp(['voxel based result has saved in output_path']);
end

%  >> After calculating and saving voxel based Myelin value in HCP or UKB cohorts:
%  >> calculate values for each brain parcels (hoacer_sn_hth parcellation file)
%--------------------------------------------------------------------------------
% -------parcellation: caluculate parcel based Myelin values------------------
%--------------------------------------------------------------------------------
disp(['parcellation has started...']);

% get mask indices
idx_parc = find(parc);
nmax = nnz(parc);
mask = parc(idx_parc);
map= map(:);

for s=1:length(subjs)

    for j = 1:max(parc(:))
        Myelin_hoacer_sn_hth(j,s)= mean(map(parc(:)==j,s), 'omitnan');
    end
end

%>>>>>>>output variable of function<<<<<<<
 Myelin_parcelvalues = Myelin_hoacer_sn_hth;
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%-------save result matrix in "path_output" ----------

save(fullfile(path_output,'Myelination_hoacer_sn_hth.mat'),'Myelin_hoacer_sn_hth','-v7.3') % Myelin matrix
disp(['parcel based result has saved in output_path']);






