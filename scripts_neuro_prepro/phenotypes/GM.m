% this script calculate Graymatter values in whole brain voxels and then
% calculate sum of GM in each parcel

%Neda SP, July 2021


% load parcellation

parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;

% 
% %Hoa parcellation
% HOA_path = '/data1/rubinov_lab/Neda/atlas';
% HOA_parc = load(fullfile(HOA_path, 'HarvardOxfordAtlas.mat')).parcellation;
% HOA_label = load(fullfile(HOA_path, 'HarvardOxfordAtlas.mat')).region_names;

%load T1w volume (bias corrected)
%>> HCP
path_input = '/data1/datasets/hcp';

%>> ukb
path_input = fullfile('/data1/datasets/ukbiobank/data/T1');
GM_path = fullfile('/T1/T1_fast');


% get file info of valid subject scans 
%>> hcp
ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5';%891-1
subjs = h5read(ts_order, '/subjects');

%>> ukb
subj_dir = dir(fullfile(path_input));
names = string({subj_dir.name});
subjs = names(3:end-2)';
%save('subjs_order.mat','subjs', '-v7.3')

%% loop through subjects
GM_vol_allvoxels = nan(length(parc(:)), length(subjs));
  N=0; 
for s = 37:length(subjs)

    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);


        % read data 

%%%>>       scan_path = fullfile(path_input, subj, 'T1w');
%         scan_name = 'T1w_acpc_dc_restore_brain.nii.gz';
%         T1w_brain = niftiread(fullfile(scan_path, scan_name));
%         %----------------------------------------
%         %We already have brain extracted data, so we skip the BET step.
%         %brain extraction (BET)
%         %bet = "/data1/rubinov_lab/bin/fsl/bin/bet";
%      
%         %command = char([bet+" "+fullfile(scan_path, 'T1w_acpc_dc_restore_brain.nii.gz')+...
%         %    " "+ fullfile(T1w_path,'T1w_acpc_dc_restore_brain.nii.gz')+...
%         %   " -m"+...
%         %   " -f .3"]);
%         %[status,cmdout] = system(command);
%         %------------------------------------------
%       
%         % Segmentation using FAST in FSL  
%         fast = "/data1/rubinov_lab/bin/fsl/bin/fast";
%%%>>        scan_path_o = fullfile('/data1/rubinov_lab/datasets/hcp/data', subj);
%         
%         command = char([fast+" -o "+fullfile(scan_path_o,'T1w_acpc_dc_restore_brain')+" "+...
%             fullfile(scan_path, 'T1w_acpc_dc_restore_brain.nii.gz')]);
%         [status,cmdout] = system(command);
%          %output: pve_0:CSF, pve_1:GM, pve_2:WM
%       
%          %loading Gray matter-->[0,1]  

        if size(dir(fullfile(path_input, subj,'T1')),1)>3
            
%          T1w_seg_GM = niftiread(fullfile(scan_path_o,'T1w_acpc_dc_restore_brain_pve_1.nii.gz'));
           T1w_seg_GM = double(niftiread(fullfile(path_input, subj, GM_path, 'T1_brain_pve_1.nii.gz')));
%         %--------
%         %resampling to MNI152_2mm
          MNI_standard_path = '/data1/rubinov_lab/bin/fsl/data/standard';
% %         MNI = niftiread(fullfile(MNI_standard_path,'MNI152_T1_2mm.nii.gz'));
            
            % first compute transform matrix using flirt in FSL
             flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";

            write_path = fullfile('/data1/rubinov_lab/Neda/GM_UKbiobank',subj);

            if ~exist(write_path)
                mkdir(write_path)
            end


            command = char([flirt+" -in "+fullfile(path_input, subj, GM_path, 'T1_brain_pve_1.nii.gz')+...
                " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                " -omat "+ fullfile(write_path,'transform_T1_to_MNI2mm.mat')]);
            [status,cmdout] = system(command);

             % apply transform
            command = char([flirt+" -in "+fullfile(path_input, subj, GM_path, 'T1_brain_pve_1.nii.gz')+...
                " -ref "+ fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
                " -omat "+ fullfile(write_path, 'transform_T1_to_MNI2mm.mat')+...
                " -out "+ fullfile(write_path, 'T1_brain_pve1_2mm.nii.gz')]);
            [status,cmdout] = system(command);


    %        GM_vol = double(niftiread(fullfile(scan_path_o,'T1w_acpc_dc_restore_brain_pve1_2mm.nii.gz')));
            GM_vol = double(niftiread(fullfile(write_path,'T1_brain_pve1_2mm.nii.gz')));      
            GM_vol = GM_vol(:);
            %for k=1:91; imshowpair(GM_vol(:,:,k),MNI(:,:,k));pause;end
            GM_vol_allvoxels(:,s) = GM_vol;

        else
            N=N+1;
            no_data(N,1) = num2cell(s);
            no_data(N,2) = cellstr(subj);

            disp([num2str(s), ': ', subj,' no data']);

        end 
end %of loop through subjs
%%
cd /data1/rubinov_lab/Neda/GM
save('GM_vol_allvoxels.mat','GM_vol_allvoxels','no_data','-v7.3') %no_data=130, all_subjs=3029

    
%% ------------------------Parcellation------------------------------------

  % parcellation
  for s=1:length(subjs)  
      
        for j = 1:max(parc(:))
            GM_115parc(j,s)= sum(GM_vol_allvoxels(parc(:)==j,s));
        end
  end
  save('GM_vol_hoacer_sn_hth.mat','GM_115parc','no_data','-v7.3')  

  
%% -------------------------ploting----------------------------------------

%----------HarvardOxford Atlas
lh_GM_1 = [GM_HOAparc(1:48,:); GM_HOAparc(97:104,:)];
    
rh_GM_1 = [GM_HOAparc(49:96,:); GM_HOAparc(105:112,:)];

figure, plot(lh_GM_1, rh_GM_1, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("HOA parc, 'Gray matter' ,  corr : "+corr(lh_GM_1(:),rh_GM_1(:))+HOA_label(1,1))
%-----115 parc
Hoacer=GM_115parc;
%lh_label indx in parc label file --> 1:48,   97:103, 114
%rh_label indx in parc label file --> 49:96, 104:110, 113

lh = [Hoacer(1:48,:); Hoacer(97:103,:); Hoacer(114,:)];
rh = [Hoacer(49:96,:);Hoacer(104:110,:); Hoacer(113,:)];

[x,y] = find(isnan(lh));
y=unique(y);

lh(:,y)=[];
rh(:,y)=[];


figure, plot(lh, rh, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("GM, UKBiobank, corr : "+corr(lh(:), rh(:)))

%% ---------121 regions parcellation 

%lh_label indx in parc label file --> 1:50,   109:116, 120
%rh_label indx in parc label file --> 51:100, 101:108, 121

lh_GMvol = [GM_vol_121parc(1:50,1);GM_vol_121parc(109:116,1); GM_vol_121parc(120,1)];
rh_GMvol = [GM_vol_121parc(51:100,1);GM_vol_121parc(101:108,1); GM_vol_121parc(121,1)];



figure, plot(nonzeros(lh_GMvol), nonzeros(rh_GMvol), '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("121parc, GM corr : "+corr(nonzeros(lh_GMvol(:), nonzeros(rh_GMvol(:)))    
    