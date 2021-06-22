%Calculating Fractional Anistrophy and mean Diffusivity using diffusion-weighted images
%first run bash file "difusion_FA_MD.sh" to generate dti__FA.nii and
%dti__MD.nii files for each subject, and then use this matlab script to
%calculate mean FA&MD in 121 brain regions

%Neda SP, June 2021


%------------------loading data--------------------------------

dataset_path = fullfile('/data1/rubinov_lab/Neda/data');
Diffusion_path = fullfile('/T1w/Diffusion');

MNI_standard_path = fullfile('/data1/rubinov_lab/bin/fsl/data/standard');
mni_2mm = uint8(niftiread(fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')));

%get a list of 1206 subjects and the paths
Files=dir(dataset_path);
Name = {Files.name};
subj_path = fullfile(dataset_path,Name)';

%% ---------------------parcellation data------------------------------------
% 
% parc_path = fullfile('/data1/rubinov_lab/Neda');
% parc_2mm = niftiread(fullfile(parc_path,'parcellation_mni.nii'));
parc_121_2mm_path = fullfile('/data1/rubinov_lab/Neda');
parc_121_2mm = niftiread(fullfile(parc_121_2mm_path, 'parcellation_mni.nii'));
parc_121_2mm(~parc_121_2mm) = nan;

%output data definition 
mean_FA = nan(max(parc_121_2mm(:)),size(subj_path,1)-2); %121*1206
mean_MD = nan(max(parc_121_2mm(:)),size(subj_path,1)-2);

%--------------------------------------------------------------------------

%load FSL
[status,cmdout] = system('module load GCC/5.4.0-2.26 OpenMPI/1.10.3 FSL/5.0.10');
N=0; %number of non-data subjects
M=0; %~^
%%
for i=609:size(subj_path,1)-2   
    %Loop through subjects
    
    temp_path = subj_path(i+2,:);
    a = numel(dir(char(fullfile(temp_path,Diffusion_path)))); %ignore directory of subjects without diffusion data
    
    if(a > 2)
        
        data = niftiread(char(fullfile(temp_path,Diffusion_path,'data.nii.gz')));
        temp_path = char(fullfile(temp_path,Diffusion_path));
        
        
        % first compute transform matrix using flirt in FSL
        flirt = "/data1/rubinov_lab/bin/fsl/bin/flirt";
        command = char([flirt+" -in "+fullfile(temp_path, 'data.nii.gz')+" -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+" -omat "+ fullfile(temp_path, 'transform.mat')]);
        [status,cmdout] = system(command);
        
        %%----------------read FA, MD and resample to 2mm------------------------
        
        FA = niftiread(string(fullfile(temp_path,'dti__FA.nii.gz')));
        MD = niftiread(string(fullfile(temp_path,'dti__MD.nii.gz')));
        
        %apply the same transform matrix to FA and MD
        command = char([flirt+" -in "+fullfile(temp_path, 'dti__FA.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -out "+fullfile(temp_path, 'dti__FA_2mm.nii.gz')+...
            " -init "+fullfile(temp_path, 'transform.mat')+...
            " -applyxfm"]);
        
        [status,cmdout] = system(command); disp(cmdout)
        
        command = char([flirt+" -in "+fullfile(temp_path, 'dti__MD.nii.gz')+...
            " -ref "+fullfile(MNI_standard_path,'MNI152_T1_2mm_brain.nii.gz')+...
            " -out "+fullfile(temp_path, 'dti__MD_2mm.nii.gz')+...
            " -init "+fullfile(temp_path, 'transform.mat')+...
            " -applyxfm"]);
        
        [status,cmdout] = system(command); disp(cmdout)
        
        disp(strcat('Resampled FA and MD for  ',char(temp_path),'  created!'))
        
        FA_2mm = niftiread(string(fullfile(temp_path,'dti__FA_2mm.nii.gz')));
        MD_2mm = niftiread(string(fullfile(temp_path,'dti__MD_2mm.nii.gz')));
        
        
        %%----calculate FA mean for each region using 121regions_2mm_parcellation-----
        
        
        for j = 1:max(parc_121_2mm(:))
            
            ix = find(parc_121_2mm==j);
            mean_FA(j,i) = mean(FA_2mm(ix));
            mean_MD(j,i) = mean(MD_2mm(ix));
               
        end
        
    else
        N=N+1;
        no_data(N,:) = cellstr(temp_path);
        disp(['Diffusion data is not existed for: ',temp_path])
        
    end
    
    %output: Fractional Anistropy(FA) & Mean Diffusivity(MD) --> [121 region * 1206 subjs]
    M=M+1;
    y_data(M,:) = cellstr(temp_path);
    X = [num2str(i),'th subject finished.'];
    disp(X)

end

%no_data=143 ---> number of subjects without diffusion data
cd /data1/rubinov_lab/Neda/
save('Diffusion_FA_MD.mat','mean_FA','mean_MD','no_data','y_data','-v7.3') % the final output matrix 

