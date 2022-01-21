%Calculating Fractional Anistrophy and mean Diffusivity using diffusion-weighted images
% 1. first run bash file "difusion_FA_MD.sh" to generate dti__FA.nii, dti__MD.nii files for each subject,
% 2. run this script to calculate mean-FA & mean-MD 

%Neda SP, 
%HCP dataset, June 2021
%Updated for UKbiobank dataset, Sep 2021


%------------------loading data--------------------------------
%>>HCP project
dataset_path = fullfile('/data1/rubinov_lab/Neda/data');
Diffusion_path = fullfile('/T1w/Diffusion');
%>>UKBiobank
dataset_path = fullfile('/data1/datasets/ukbiobank/data/diffusion');
Diffusion_path = fullfile('/dMRI/dMRI');


MNI_standard_path = fullfile('/data1/rubinov_lab/bin/fsl/data/standard');
mni_2mm = uint8(niftiread(fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')));


% get file info of valid subjects
%>> hcp
ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5';%891
subjs = h5read(ts_order, '/subjects');

%>> ukb
subj_dir = dir(fullfile(dataset_path));
names = string({subj.name});
subjs = names(3:end-2)';
%save('subjs_order.mat','subjs', '-v7.3')

% load parcellation

parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;

%--------------------------------------------------------------------------

%load FSL module in Platypus
[status,cmdout] = system('module load GCC/5.4.0-2.26 OpenMPI/1.10.3 FSL/5.0.10');
N=0; %number of no data subjects

%% Loop through subjects
%N_voxels = 91*109*91;
FA_allvoxels = nan(length(parc(:)), length(subjs));
MD_allvoxels = nan(length(parc(:)), length(subjs));
        

for i=1:length(subjs)   
   
    subj_i = num2str(subjs(i));
    disp([num2str(i), ': ', subj_i]);
    
    %check for availibility of diffusion data-->just HCP
%    if size(dir(fullfile(dataset_path,subj_i,Diffusion_path)),1)>3

         %>>>>>>>>>>>>>commented lines are already done<<<<<<<<<<<<<<<<<<<<<       
               
        %load diffusion data
%        data = niftiread(char(fullfile(dataset_path, subj_i, Diffusion_path,'data.nii.gz')));

%         % first compute transform matrix using flirt in FSL
%         flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
%>> hcp
%          temp_path = fullfile(dataset_path, subj_i,Diffusion_path);
%         command = char([flirt+" -in "+fullfile(temp_path, 'data.nii.gz')+...
%             " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%             " -omat "+ fullfile(temp_path, 'transform.mat')]);
%         [status,cmdout] = system(command);

%>> ukb
    if size(dir(fullfile(dataset_path,subj_i,Diffusion_path)),1)>3
        
        temp_path = fullfile(dataset_path, subj_i, Diffusion_path);
        write_path = fullfile('/data1/rubinov_lab/Neda/UKbiobank',subj_i);
        
        if ~exist(write_path)
            mkdir(write_path)
        end
        
        
        %FA
        command = char([flirt+" -in "+fullfile(temp_path, 'dti_FA.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -omat "+fullfile(write_path, 'dti_FA_to_MNI2mm.mat')]);
        [status,cmdout] = system(command);
        
        %MD
         command = char([flirt+" -in "+fullfile(temp_path, 'dti_MD.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -omat "+fullfile(write_path, 'dti_MD_to_MNI2mm.mat')]);
        [status,cmdout] = system(command);

        %%----------------read FA, MD (outout of "difusion_FA_MD.sh" script) and resample to 2mm------------------------
%>> hcp        
%         FA = niftiread(string(fullfile(temp_path,'dti__FA.nii.gz')));
%         MD = niftiread(string(fullfile(temp_path,'dti__MD.nii.gz')));
%>> ukb
%         %FA = niftiread(char(fullfile(dataset_path, Diffusion_path, subj_i,'/dMRI/dMRI/dti_FA.nii.gz')));
%         FA_MNI_1mm = niftiread(char(fullfile(dataset_path, Diffusion_path, subj_i,'/dMRI/TBSS/FA/dti_FA_to_MNI.nii.gz')));
%         % load data and resample to MNI_2mm
%         [x, y, z] = size(FA_MNI_1mm);
%         [Xq, Yq, Zq] = ndgrid(1.5:2:x, 1.5:2:y, 1.5:2:z);
%         FA_MNI_2mm = double(interpn(FA_MNI_1mm, Xq, Yq, Zq));

%>>>>>>>>>apply the same transform matrix to FA and MD

%>> hcp        
%         command = char([flirt+" -in "+fullfile(temp_path, 'dti__FA.nii.gz')+...
%             " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
%             " -out "+fullfile(temp_path, 'dti__FA_2mm.nii.gz')+...
%             " -init "+fullfile(temp_path, 'transform.mat')+...
%             " -applyxfm"]);

%>>ukb
        %FA
        command = char([flirt+" -in "+fullfile(temp_path, 'dti_FA.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -out "+fullfile(write_path, 'dti_FA_2mm.nii.gz')+...
            " -init "+fullfile(write_path, 'dti_FA_to_MNI2mm.mat')+...
            " -applyxfm"]);
        [status,cmdout] = system(command); disp(cmdout)
        
        %MD
        command = char([flirt+" -in "+fullfile(temp_path, 'dti_MD.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -out "+fullfile(write_path, 'dti_MD_2mm.nii.gz')+...
            " -init "+fullfile(write_path, 'dti_MD_to_MNI2mm.mat')+...
            " -applyxfm"]);
        [status,cmdout] = system(command); disp(cmdout)
        
        disp(strcat('Resampled FA and MD for   ', char(subj_i),'  created!'))
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>> hcp
%        %load resampled FA, MD data
%         FA_2mm = niftiread(string(fullfile(temp_path,'dti__FA_2mm.nii.gz')));
%         MD_2mm = niftiread(string(fullfile(temp_path,'dti__MD_2mm.nii.gz')));
%          %output: Fractional Anistropy(FA) & Mean Diffusivity(MD)
%         FA_allvoxels(:,i) = FA_2mm(:);
%         MD_allvoxels(:,i) = MD_2mm(:);
%         disp([num2str(i), ': ', subj_i, ' finished']);
%>> ukb
        %load resampled FA, MD data
        FA_2mm = double(niftiread(string(fullfile(write_path,'dti_FA_2mm.nii.gz'))));
        MD_2mm = double(niftiread(string(fullfile(write_path,'dti_MD_2mm.nii.gz'))));
         %output: Fractional Anistropy(FA) & Mean Diffusivity(MD)
         %for k=1:91; imshowpair(MD_2mm(:,:,k),mni_2mm(:,:,k));pause;end
        FA_allvoxels(:,i) = FA_2mm(:);
        MD_allvoxels(:,i) = MD_2mm(:);
        disp([num2str(i), ': ', subj_i, ' finished']);
%         
     else
        N=N+1;
        no_data(N,1) = num2cell(i);
        no_data(N,2) = cellstr(subj_i);
        
        disp([num2str(i), ': ', subj_i,' no data']);
     end   
end


function myniftiread(filename)

%% save 
%no_data=143 ---> number of subjects without diffusion data
cd /data1/rubinov_lab/Neda/Diffusion
save('Diffusion_FA_MD_allvoxels.mat','FA_allvoxels','MD_allvoxels','no_data','-v7.3') % the final output matrix 
%>> ukb
cd /data1/rubinov_lab/Neda/Diffusion_UKB
save('Diffusion_FA_MD_allvoxels.mat','FA_allvoxels','MD_allvoxels','-v7.3')
%% ------------------------Parcellation------------------------------------
%Hoa parcellation
HOA_path = '/data1/rubinov_lab/Neda/atlas';
HOA_parc = load(fullfile(HOA_path, 'HarvardOxfordAtlas.mat')).parcellation;
HOA_label = load(fullfile(HOA_path, 'HarvardOxfordAtlas.mat')).region_names;

  %% parcellation
  for s=1:length(subjs)  
      
        for j = 1:max(parc(:))
            FA_115parc(j,s)= nanmean(FA_all(parc(:)==j,s));
            MD_115parc(j,s)= nanmean(MD_all(parc(:)==j,s));
        end
  end
  
 save('Diffusion_FA_MD_hoacer_sn_hth.mat','FA_115parc','MD_115parc','no_data','-v7.3') % the final output matrix  
  
%% -------------------------ploting----------------------------------------

%----------HarvardOxford Atlas
%---FA---
lh_GM_1 = [FA_HOAparc(1:48,:); FA_HOAparc(97:104,:)];
    
rh_GM_1 = [FA_HOAparc(49:96,:); FA_HOAparc(105:112,:)];

%remove nans from matrix
[~,y] = find(isnan(lh_GM_1));y=unique(y);
lh_GM_1(:,y)=[];
rh_GM_1(:,y)=[];

figure, plot(lh_GM_1, rh_GM_1, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("HOA parc, 'mean Fractional Anisotrophy' ,  corr : "+corr(lh_GM_1(:),rh_GM_1(:)))


%---MD---
lh_GM_1 = [MD_HOAparc(1:48,:); MD_HOAparc(97:104,:)];
    
rh_GM_1 = [MD_HOAparc(49:96,:); MD_HOAparc(105:112,:)];

%remove nans from matrix
[x,y] = find(isnan(lh_GM_1));y=unique(y);
lh_GM_1(:,y)=[];
rh_GM_1(:,y)=[];

figure, plot(lh_GM_1, rh_GM_1, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("HOA parc, 'mean Mean Diffusivity' ,  corr : "+corr(lh_GM_1(:),rh_GM_1(:)))


%------  115
%lh_label indx in parc label file --> 1:48,   97:103, 114
%rh_label indx in parc label file --> 49:96, 104:110, 113

hoacer = MD_115parc;

lh_hoacer = [hoacer(1:48,:); hoacer(97:103,:); hoacer(114,:)];
rh_hoacer = [hoacer(49:96,:);hoacer(104:110,:); hoacer(113,:)];

% remove nans
[x1,y1] = find(isnan(lh_hoacer)); 
y1=unique(y1);
lh_hoacer(:,y1)=[];
rh_hoacer(:,y1)=[];

% for k=1:length(y1)
%     number=y1(k,1);
%     no = subjs(number,1);
%     y1(k,2) = no;
% end


figure, plot(lh_hoacer, rh_hoacer, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("MD, UKBiobank, lh/rh corr : "+corr(lh_hoacer(:), rh_hoacer(:)))
        