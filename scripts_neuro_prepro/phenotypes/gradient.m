%Running on each subject separately

% Diffusion embedding scrpit for HCP dataset
% W: preproccesed clean timeseries for all voxels
% This script call >> diffusion_embedding.m function located at : "cd /data1/rubinov_lab/Neda/mika"

timeseries_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_voxels/timeseries';
% subj_fMRI_dir = dir(fullfile(timeseries_path));
% names = string({subj_fMRI_dir.name})';
% subjs = names(3:end-1);

% get file info of valid subject timeseries 
ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5';%890
subjs = h5read(ts_order, '/subjects');


%loading mask for nonzero voxels
mask = load(fullfile('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_voxels/','mask.mat')).mask(:);
mask_ind = find(mask);

% number of scans in HCP
scans=4;

% load parcellation
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;
lh_lables = [1:48];
rh_lables= [49:96];
subcort_lables= [97:115];

%get index for lh, rh, subcortical regions

idx_left = false(size(parc));
for i=1:48
    idx_left = idx_left | (parc == i);
end

idx_right = false(size(parc));
for i=49:96
    idx_right = idx_right | (parc == i);
end

idx_subcort = false(size(parc));
for i=97:115
    idx_subcort = idx_subcort | (parc == i);
end

idx_left = idx_left(mask);
idx_right = idx_right(mask);
idx_subcort = idx_subcort(mask);
%%

 U1_brain=cell(length(subjs),3);

for s = 1:length(subjs)       %through subjects
    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    
    
    %load data : mat files
    mat_name = [subj '.mat'];   
    temp_path = fullfile(timeseries_path, mat_name);
    data = load(temp_path);
    
    t= table2array(data.scans(1,2)); 
    %W_scans=nan(length(mask_ind),t,scans);
    
    next=1;
    
    lh_scans = zeros(nnz(idx_left),nnz(idx_left));  lh_freq_scans=lh_scans; 
    rh_scans = zeros(nnz(idx_right),nnz(idx_right)); rh_freq_scans=rh_scans;
    subcort_scans = zeros(nnz(idx_subcort),nnz(idx_subcort)); sub_freq_scans=subcort_scans;

    for ii=1:scans       % loop through fMRI scan metrices
        timeseries = cell2mat(data.Vp_clean(ii,1));

        
        if  isempty(timeseries)==1 %SKIP empty scans
            X = [num2str(s),': ',subj,' ',num2str(ii),'th scan is empty'];disp(X)
            no_scan(next,:) = cellstr(X);
            next=next+1;
            continue
            
        else
            
             %calculate correlation matrix for each scan        
             lh_corr = corr(timeseries(idx_left, :)');
             rh_corr = corr(timeseries(idx_right,:)');
             subcort_corr = corr(timeseries(idx_subcort,:)');
             %all scans
             lh_scans = lh_scans + lh_corr;  %sum of correlation matrices
             rh_scans = rh_scans + rh_corr;
             subcort_scans = subcort_scans + subcort_corr;
             %# of repeats for non-infinite values 
             lh_freq_scans = lh_freq_scans + isfinite(lh_corr); 
             rh_freq_scans = rh_freq_scans + isfinite(rh_corr);
             sub_freq_scans = sub_freq_scans + isfinite(subcort_corr);
        end 
     
    end %scans   
    
    % compute scans average here for each subjects
    lh_corr_avg = lh_scans ./ lh_freq_scans;
    rh_corr_avg = rh_scans ./ rh_freq_scans;
    subcort_corr_avg = subcort_scans ./ sub_freq_scans;
      

 %%  pass data to function
 
    data_size = {size(lh_corr_avg,1), size(rh_corr_avg,1), size(subcort_corr_avg,1)};
    
   % Get indices of finite values (non NaNs)
   [xf_l,~]=find(isfinite(lh_corr_avg(:,1:20)));
   [xf_r,~]=find(isfinite(rh_corr_avg(:,1:20)));
   [xf_s,~]=find(isfinite(subcort_corr_avg(:,1:20)));
   
   ix_finite = {unique(xf_l), unique(xf_r), unique(xf_s)};
   
   % Get indices of infinite values (NaNs)
   
   xi_l = setdiff([1:size(lh_corr_avg,1)]',xf_l);
   xi_r = setdiff([1:size(rh_corr_avg,1)]',xf_r);
   xi_s = setdiff([1:size(subcort_corr_avg,1)]',xf_s);
  
   %exclude nans --> prepare function input without nans
  subcort_corr_avg(unique(xi_s),:)=[]; subcort_corr_avg(:,unique(xi_s))=[];  
  lh_corr_avg(unique(xi_l),:)=[]; lh_corr_avg(:,unique(xi_l))=[];
  rh_corr_avg(unique(xi_r),:)=[]; rh_corr_avg(:,unique(xi_r))=[];
  

    %-----------------compute diffusion embedding ------------------------
    % inputs: W= timeseries matrix, c= number of output gradients, r= number of computed singular vectors (to speed up)
    % output: U1= gradient maps
    
    disp([subj, 'gradient map calculation started ...']);
    c=1;
    W_brain = {lh_corr_avg, rh_corr_avg, subcort_corr_avg}; % contains voxel timesries of lh, rh, and subcortical separately
    
    for j=1:length(W_brain)      
        U1 = nan(data_size{1,j},1);
        U1(ix_finite{1,j}) = diffusion_embedding(W_brain{1,j}, c);
        
        U1_brain{s,j} = U1; % a 1x3 cell for each subject
    end
    
    disp([num2str(s), ' : ', subj ' finished']);
    
end %subjects

cd /data1/rubinov_lab/Neda/Gradient
save('Gradient_allsubjects.mat','U1_brain','-v7.3')
disp('result has been saved!');
%result have been updated DEC 10
%% LOAD DATA

U1_brain = load('Gradient_allsubjects.mat');
U1_brain_avg = load('Gradient_allsubjs_avg.mat').U1_brain;
    
    %reshape the data
    vol_avg =zeros(91,109,91);
    vol_avg(idx_left) = U1_brain_avg{1,1};
    vol_avg(idx_right) = U1_brain_avg{1,2};
    vol_avg(idx_subcort) = U1_brain_avg{1,3};

%% %%
%>> visualize and parcellate each subject
for s = 1:length(subjs)       %through subjects
    
    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    
    %reshape the data
    vol =zeros(91,109,91);
    vol(idx_left) = U1_brain{s,1};
    vol(idx_right) = U1_brain{s,2};
    vol(idx_subcort) = U1_brain{s,3};
    vol = vol(:);
    
    %for i = 1:91; imagesc(vol(:, :, i)); axis equal; axis off; pause; end
    
    for j = 1:max(parc(:))
        Gradient_115parc(j,s)= nanmean(nonzeros(vol(parc(:)==j)));
    end
        
end
cd /data1/rubinov_lab/Neda/Gradient
save('HCPgradient_hoacer_sn_hth.mat','Gradient_115parc','-v7.3') %[115,890]
disp('result has been saved!');

%% load allsubj data in volume format
%Transfer each individual into allsubjs_average gradient map.
for s = 1:length(subjs)       %through subjects
    
    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    
    %reshape the data
    vol =zeros(91,109,91);
    vol(idx_left) = U1_brain{s,1};
    vol(idx_right) = U1_brain{s,2};
    vol(idx_subcort) = U1_brain{s,3};
    
    %svd for each brain slice
    for i = 1:91
        %imagesc(vol(:, :, i));axis equal; axis off; pause;
        
        temp = vol_avg(:,:,i).*vol(:,:,i);
        %replace nan with zero
        idx_n = find(isnan(temp));temp(idx_n)=0;
        
        [S,~,T] = svd (temp,'econ');
        vol(:,:,i) = vol(:,:,i).*(T*S')';
               
    end
    for i = 1:91; imagesc(vol(:, :, i)); axis equal; axis off; pause; end

 
        
end



%% plot
Gradient_Hoacer=Gradient_115parc;
%lh_label indx in parc label file --> 1:48,   97:103, 114
%rh_label indx in parc label file --> 49:96, 104:110, 113

lh_Gradient = [Gradient_Hoacer(1:48,:); Gradient_Hoacer(97:103,:); Gradient_Hoacer(114,:)];
rh_Gradient = [Gradient_Hoacer(49:96,:);Gradient_Hoacer(104:110,:); Gradient_Hoacer(113,:)];

[x,y] = find(isnan(rh_Gradient));
y=unique(y);

lh_Gradient(:,y)=[];
rh_Gradient(:,y)=[];


figure, plot(lh_Gradient, rh_Gradient, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("Gradient, HCP, corr : "+corr(lh_Gradient(:), rh_Gradient(:)))

