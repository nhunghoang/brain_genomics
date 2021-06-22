%This script calculate the mean of Myelination in specific brain regions 
%Here we use a parcellation with 121 brain regions

%Neda SP, June 2021

%***note: this script should be modified for all HCP subjects


%---------------------Loading data-----------------------------------------
path = '/data/rubinov_lab/datasets/hcp/data/single_subject_full/100307/T1w';
Myelination = niftiread(fullfile(path,'T1wDividedByT2w.nii.gz'));

MNI_standard_path = fullfile('/data1/rubinov_lab/bin/fsl/data/standard');
mni_2mm = uint8(niftiread(fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')));


%% ---------------------parcellation data------------------------------------
% 
parc_121_2mm_path = fullfile('/data1/rubinov_lab/Neda');
parc_121_2mm = niftiread(fullfile(parc_121_2mm_path, 'parcellation_mni.nii'));
parc_121_2mm(~parc_121_2mm) = nan;

%output data definition 
%mean_Myelin = nan(max(parc_121_2mm(:)),size(subj_path,1)-2); %121*1206


%% -----------resampling data to 2mm MNI space using FSL command----------------

% first compute transform matrix using flirt in FSL
        flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
        command = char([flirt+" -in "+fullfile(path, 'T1wDividedByT2w.nii.gz')+" -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+" -omat "+ fullfile(path, 'transform_to_MNI2mm.mat')]);
        [status,cmdout] = system(command);
        
%apply the same transform matrix to input
        command = char([flirt+" -in "+fullfile(path, 'T1wDividedByT2w.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -out "+fullfile(path, 'T1wDividedByT2w_MNI2mm.nii.gz')+...
            " -init "+fullfile(path, 'transform_to_MNI2mm.mat')+...
            " -applyxfm"]);
        
        [status,cmdout] = system(command); disp(cmdout)
                
%loading transformed file
Myelination_MNI2mm = niftiread(fullfile(path,'T1wDividedByT2w_MNI2mm.nii.gz'));



%%----calculate Myelination mean for each region using 121regions_2mm_parcellation-----

%for i=1:size(subj_path,1)-2  
    i=1;
        for j = 1:max(parc_121_2mm(:))
            
            ix = find(parc_121_2mm==j);
            mean_Myelin(j,i) = mean(Myelination_MNI2mm(ix));
               
        end
        
    %output: Myelination value for all 121 region --> [121 region * 1 subj]

%end


                