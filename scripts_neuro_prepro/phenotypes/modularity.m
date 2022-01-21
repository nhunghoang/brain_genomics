% this script calculate the modularity of brain network 
% modified for HCP dataset
%last update January 2022


%% calculate corr matrices of all subjects
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


corr_allsubjs=cell(length(subjs),1);

for s = 1:length(subjs)       %through subjects
    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    
    corr_scans=zeros(max(parc),max(parc));freq_scans=corr_scans;
    
    %load data : mat files
    mat_name = [subj '.mat'];   
    temp_path = fullfile(timeseries_path, mat_name);
    data = load(temp_path);
    
    t= table2array(data.scans(1,2)); 
    %W_scans=nan(length(mask_ind),t,scans);
    
    next=1;
    for ii=1:scans       % loop through fMRI scan metrices
        timeseries = cell2mat(data.Vp_clean(ii,1));
  
        if  isempty(timeseries)==1 %SKIP empty scans
            X = [num2str(s),': ',subj,' ',num2str(ii),'th scan is empty'];disp(X)
            no_scan(next,:) = cellstr(X);
            next=next+1;
            continue
            
        else
            masked_parc=parc(mask);
            %parcellation
            for j = 1:max(parc(:))
                temp=timeseries(masked_parc==j,:);
                [x_0,~]=find(temp==0);
                temp(unique(x_0),:)=[];
                V_115parc(j,:) = nanmean(temp);
            end 
            
             %calculate correlation matrix for each scan        
             corr_mat = corr(V_115parc');

             %all scans
             corr_scans = corr_scans + corr_mat;  %sum of correlation matrices
 
             %# of repeats for non-infinite values 
             freq_scans = freq_scans + isfinite(corr_mat); 
        end 
     
    end %scans   
    
    % compute scans average here for each subjects
    corr_avg = corr_scans ./ freq_scans;
    
    %STORE DATA of all subjects 
    corr_allsubjs{s,1}=corr_avg;

 end %subjects
 
    
%     % save average of all subjects's corr matrix
cd /data1/rubinov_lab/Neda/community
save('corr_HCP_parc115.mat','corr_allsubjs','-v7.3')
disp('result has been saved!');
    
 %%  pass data to function
 %load corr matrices
 corr_allsubjs = load('corr_HCP_parc115.mat').corr_allsubjs;
 
cell_W=corr_allsubjs;

gamm_range = 1:0.02:2;

%% remove the main diagonal in al matrices

for i = 1:numel(cell_W)
    n=size(cell_W{i},2);
    cell_W{i}(1:n+1:end) = 0;
end

%%

M = multiscale_community(cell_W, gamm_range);

%number of modules detected in each run
figure, plot(max(M, [], 1))
xlabel('gamma'); ylabel('number of communities');axis square

figure
[~,NMI] = partition_distance(M);
imagesc(NMI>0.99); axis square, colorbar
xlabel('partition'), ylabel('partition')
title('NMI>0.99')
%imagesc(NMI(10:15,10:15)>0.99)

%% pairs of nodes in the same module (within-module probability) 
% 0 : pairs of nodes in different module,  1 : pairs of nodes in the same module
% convert the module affiliation vector into a binary matrix

P = zeros (115);
%for i = 1:numel(cell_W)
    
    for j=1:size(M,2)       
        P = P + (M(:,j)==M(:,j)');
    end
    %probability
    P = P/size(M,2);
    imagesc(P); axis square; colorbar
    xlabel('node')
    ylabel('node')
    title('within-module probability')
    
%end

plot(max(M, [], 1))
tabulate(M(:, 15))
tabulate(M(:, 20))
addpath /home/rubinom/BCT/2017_01_15_BCT/
% This function reorders the connectivity matrix by modular structure and may consequently be useful in visualization of modular structure.
ix = reorder_mod(P, M(:, 20));



imagesc(P(ix,ix))
xlabel(parc_label(ix,:))

%% calculate the correlation of nodes within modules

    
corr_module_allsubjects = cell(length(subjs),max(M(:,20)));

for s = 1:length(subjs)       %through subjects
    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
   
    
    %load data : mat files
    mat_name = [subj '.mat'];   
    temp_path = fullfile(timeseries_path, mat_name);
    data = load(temp_path);
    
    t= table2array(data.scans(1,2)); 
    %W_scans=nan(length(mask_ind),t,scans);
    
    
    V_scans = zeros(max(parc),1200); V_freq_scans=V_scans;
    
    next=1;
    for ii=1:scans       % loop through fMRI scan metrices
        timeseries = cell2mat(data.Vp_clean(ii,1));
  
        if  isempty(timeseries)==1 %SKIP empty scans
            X = [num2str(s),': ',subj,' ',num2str(ii),'th scan is empty'];disp(X)
            no_scan(next,:) = cellstr(X);
            next=next+1;
            continue
            
        else
            masked_parc=parc(mask);
            %parcellation
            for j = 1:max(parc(:))
                temp=timeseries(masked_parc==j,:);
                [x_0,~]=find(temp==0);
                temp(unique(x_0),:)=[];
                V_115parc(j,:) = nanmean(temp);
            end 
            
             %all scans
             V_scans = V_scans + V_115parc;  %sum of correlation matrices
             %# of repeats for non-infinite values 
             V_freq_scans = V_freq_scans + isfinite(V_scans);   
        end 
     
    end %scans   
    % compute scans average here for each subjects
    V_scans_avg = V_scans ./ V_freq_scans;
              
    %>>>>>>>>>>extract  within module nodes<<<<<<<<<<<<<<<<<<<<
     for m=1:max(M(:,20))
         idx = find(M(:, 20) == m);
         module_timeseries = V_scans_avg(idx,:);
                 
         %calculate correlation matrix for each node        
         corr_module = corr(module_timeseries');
         %STORE DATA of all subjects 
         corr_module_allsubjects{s,m} = corr_module;
                 
     end

 end %subjects
 
 %visualization
 for m=1:max(M(:,20))
    temp =cell2mat(corr_module_allsubjects(1,m));
    figure, imagesc(temp), colorbar, axis square, hold on   
    xlabel(length(temp)+" nodes"), ylabel(length(temp)+" nodes"),
    title("module : "+m)
    
    
 end
 
 