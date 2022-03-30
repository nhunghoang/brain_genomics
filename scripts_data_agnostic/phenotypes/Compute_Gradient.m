

function [Avg_Gradient, allsubjs_Gradient, allsubjs_Gradient_svd] = Compute_Gradient(project, path_output)

%-----------------------------------------------------------------
%
%   function summary:
%     this script load timeseries matrices (.mat file) as input in forms of: [parcels, timepoints] or [allvoxels,timepoints]
%       >>>>> processing steps:
%        >parcel timeseries (115 parcellation) and calculate correlation matrix for each scan
%	     >each subjects: average corr matrix of all scans
%	     >calculation group average correlation matrix in that cohort
%	     >pass correlation matrices(group avg & each subject separately) to "diffusion embedding" funstion
%	     >transfer gradient matrix of each subject to the group avg gradient matrix
%
%   function inputs:
%     project:      label: human connectome project("hcp") or UK-biobank("ukb")
%     path_output:  path to store output Like: path_output = '/data1/rubinov_lab/Neda/Gradient_UKB'
%
%   function outouts:
%     Avg_Gradient        :  group average gradient map :[115, 1]
%     allsubjs_Gradient   :  gradient map of all subjects separately: [115, #subjects]
%     allsubjs_Gradient_svd : gradient map of all subjects separately, transfered into group avg map: [115, #subjects]
%
%  The result matrices for "HCP" and "UKB" are saved in the following paths:
%---------UKB------------- 
%  path_output = '/data1/rubinov_lab/Neda/Gradient_UKB/'
%   >> GradientAvg_115corr.mat; 
%   >> Gradient_allsubjs.mat;
%   >> Gradient_allsubjs_to_avg.mat;

%---------HCP-----------
%  path_output = '/data1/rubinov_lab/Neda/Gradient/'
%   >> Gradient_allsubjs_avg.mat; 
%   >> Gradient_HCP_hoacer_sn_hth.mat;
%   >> Gradient_HCP_hoacer_sn_hth_svd.mat;
%
%*Note: the calculation pipeline for both cohorts is the same!
%
% Neda Sardaripour, 2021 
%-----------------------------------------------------------------
%-----------------------------------------------------------------

%----- load parcellation: hoacer_sn_hth (115 brain regions) ------
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/2021_hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;
% lh_lables = [1:48];
% rh_lables= [49:96];
% subcort_lables= [97:115];

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

        
%loading mask for nonzero voxels
mask = nonzeros(parc);
mask_ind = find(mask);

idx_left = idx_left(mask);
idx_right = idx_right(mask);
idx_subcort = idx_subcort(mask);



%------------------------------------------------------------------------------------------       
%>>>>> switch between two projects: 1. Human Connectome Project(hcp) 2.UKBiobank(ukb) <<<<<%
%------------------------------------------------------------------------------------------
switch project
    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>Human Connectome Project (HCP)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------
    case "hcp"

        disp(['Gradient calculation in HCP dataset has started...']);
        % input path
        timeseries_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_voxels/timeseries';

        % get file info of valid subject scans
        ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5';%890
        subjs = h5read(ts_order, '/subjects');

        %output variables
        %storing average of all subject corr matrix
        corr_all=zeros(max(parc(:)),max(parc(:)));  freq_all=corr_all;
        next=1; % for storing no data subjects

        % storing all subjects corr matrix separately
        U1_brain=cell(length(subjs),1);

        %loop through subjects
        for s = 1:length(subjs)
            subj = int2str(subjs(s));
            disp([num2str(s), ': ', subj]);


            %load data : mat files
            mat_name = [subj '.mat'];
            temp_path = fullfile(timeseries_path, mat_name);
            data = load(temp_path);


            all_scans = zeros(max(parc(:)),max(parc(:)));  all_freq_scans=all_scans;
            timeseries_115parc = zeros(max(parc(:)), 1200);

            for ii=1:4       % loop through fMRI scan metrices
                timeseries = cell2mat(data.Vp_clean(ii,1));

                if  isempty(timeseries)==1 %SKIP empty scans
                    X = [num2str(s),': ',subj,' ',num2str(ii),'th scan is empty'];disp(X)
                    no_scan(next,:) = cellstr(X);
                    next=next+1;
                    continue
                else

                    % -------parcellation------------------
                    for j = 1:max(parc(:))
                        temp = timeseries(mask==j,:);
                        [x,~] =  find(temp(:,1)== 0);
                        temp(unique(x),:)=[];
                        timeseries_115parc(j,:) = mean(temp,'omitnan');
                    end
                end

                %calculate correlation matrix for each parcellated scan
                scan_corr = corr(timeseries_115parc');
                %all scans
                all_scans = all_scans + scan_corr;  %sum of correlation matrices

                %# of repeats for non-infinite values
                all_freq_scans = all_freq_scans + isfinite(scan_corr);

            end %scans

            % compute scans average here for each subjects
            corr_subj = all_scans ./ all_freq_scans;
            % storing data for all subjects
            U1_brain{s,:}=corr_subj;

            %.........this part added for calculating avg of all subject..........%
            %STORE DATA of all subjects
            %averaging
            corr_subj(isnan(corr_subj))=0;
            corr_all = corr_all + corr_subj;
            freq_all = freq_all + isfinite(corr_subj);
            %......................................................................%
            mm=0;
        end %subjects

        % compute corr average of all subjects here
        corr_avgall = corr_all ./ freq_all;
        

        %-------save correlation matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/Gradient
        save(fullfile(path_output,'GroupAvg_115corr.mat'),'corr_avgall','-v7.3')
        save(fullfile(path_output,'corr_allsubjs_hoacer_sn_hth.mat.mat'),'U1_brain','-v7.3')
        display('correlation matrices has been saved');

        %--------------------------------------------------------------------------------------
        %       >>>>>>>>>>>>>>>>>>>>UK Biobank (UKB)<<<<<<<<<<<<<<<<<<<<<
        %--------------------------------------------------------------------------------------
        %Feb 15, 2021 updated With new final_subjs.txt order file

    case "ukb"

        disp(['Gradient calculation in ukb dataset has started...']);
        % input path
        path_input = '/data1/rubinov_lab/brain_genomics/data_UKB/hoacer_sn_hth_voxels/timeseries';
        timeseries_path = path_input;
        %path_output = '/data1/rubinov_lab/Neda/Gradient_UKB/'


        % get list of subjects
        UKB_subjs_order = fullfile('/data1/rubinov_lab/brain_genomics/data_UKB');
        subjs = load(fullfile(UKB_subjs_order, 'ordered_subject_list.txt'));


        %output variables
        %storing average of all subject corr matrix
        corr_all=zeros(max(parc(:)),max(parc(:)));  freq_all=corr_all;
        next=1;
        % storing all subjects corr matrix separately
        U1_brain=cell(length(subjs),1);

        %loop through subjects
        for s = 1:length(subjs)
            subj = int2str(subjs(s));
            disp([num2str(s), ': ', subj]);


            %load data : mat files
            mat_name = [subj];
            temp_path = fullfile(timeseries_path, mat_name);
            data = load(temp_path);


            all_scans = zeros(max(parc(:)),max(parc(:)));  all_freq_scans=all_scans;
            timeseries_115parc = zeros(max(parc(:)), data.scans.tmax);

            scan=1;
            for ii=1:scan
                timeseries = cell2mat(data.Vp_clean(ii,1));

                if  isempty(timeseries)==1 %SKIP empty scans
                    X = [num2str(s),': ',subj,' ',num2str(ii),'th scan is empty'];disp(X)
                    no_scan(next,:) = cellstr(X);
                    next=next+1;
                    continue
                else

                    % -------parcellation------------------

                    for j = 1:max(parc)
                        temp = timeseries(mask==j,:);
                        [x,~] =  find(temp(:,1)== 0);
                        temp(unique(x),:)=[];
                        timeseries_115parc(j,:) = mean(temp, 'omitnan');
                    end
                end

                %calculate correlation matrix for each parcellated scan
                scan_corr = corr(timeseries_115parc');
                %all scans
                all_scans = all_scans + scan_corr;  %sum of correlation matrices

                %# of repeats for non-infinite values
                all_freq_scans = all_freq_scans + isfinite(scan_corr);

            end %scans

            % compute scans average here for each subjects
            corr_subj = all_scans ./ all_freq_scans;
            % storing data for all subjects
            U1_brain{s,:}=corr_subj;

            %.........this part added for calculating avg of all subject..........%
            %STORE DATA of all subjects
            %averaging
            corr_subj(isnan(corr_subj))=0;
            corr_all = corr_all + corr_subj;
            freq_all = freq_all + isfinite(corr_subj);
            %......................................................................%

        end %subjects

        % compute corr average of all subjects here
        corr_avgall = corr_all ./ freq_all;


        %-------save correlation matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/Gradient
        save(fullfile(path_output,'GroupAvg_115corr.mat'),'corr_avgall','-v7.3')
        save(fullfile(path_output,'corr_allsubjs_hoacer_sn_hth.mat'),'U1_brain','-v7.3')
        display('correlation matrices has been saved');
end


%------------------------------------------------------------------------%
%----------------- Gradient, (diffusion embedding algorithm) ------------%
%------------------------------------------------------------------------%
%  >> After calculating and saving correlation matrices

%U1_brain -> contains corr matrix for all subjects***********change it
U1_brain=load(fullfile(path_output,'corr_allsubjs_hoacer_sn_hth.mat')).U1_brain;
%corr_avgall -> contains one group average corr matrix 
corr_avgall=load(fullfile(path_output,'GroupAvg_115corr.mat')).corr_avgall;

data = {{corr_avgall}; U1_brain}; %2*1 cell

%% run over two sets of data: group avg & all subjects separately. 
for j=1:size(data,1)

    data_temp = data{j};

    for i =1:size(data_temp,1)

        data1 = data_temp{i};

        % fullzeros rows and column are replced with  Nan values (because these values were originally NaN
        % : nan values were replaced with zero when correlation matrix was calculated in "line 256")
        % this is because in parcellation some parcells do not have
        % timeseries (limited number of scans from each subject)
        zero_all=(data1(:,1)==0);
        data1(zero_all,:)=NaN;data1(:,zero_all)=NaN;

        data_size = size(data1,1);

        % Get indices of finite values (non NaNs)
        [xf_l,~]=find(isfinite(data1(:,1:20)));
        ix_finite = {unique(xf_l)};

        % Get indices of infinite values (NaNs)
        xi_l = setdiff([1:size(data1,1)]',xf_l);

        %exclude nans --> prepare function input without nans
        data1(unique(xi_l),:)=[]; data1(:,unique(xi_l))=[];

        %-----------------pass to function ------------------------
        % follows the Marguiles et al. (2016) approach.

        % inputs: W= timeseries matrix, c= number of output gradients, r= number of computed singular vectors (to speed up)
        % output: U1= gradient maps

        disp([num2str(i), ':  , 1st gradient ...']);
        c=1;
        W_brain = {data1};
        %-----------------------------------------------------------
        for jj=1:length(W_brain)

            U1 = nan(data_size,1);
            U1(ix_finite{1,jj}) = diffusion_embedding(W_brain{1,jj}, c);

            %storing result
            %>>>>>>>output variable of function<<<<<<<
            if j==1
                Avg_Gradient = U1;
            else
                %all subjs :{890,1}
                allsubjs_Gradient{i,jj} = U1;
            end
            %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        end
        %-----------------------------------------------------------
    end

    disp([num2str(j), ' :  finished']);

end

%save data in output path
save(fullfile(path_output,'Gradient_GroupAvg.mat'),'Avg_Gradient','-v7.3')
save(fullfile(path_output,'Gradient_allsubjs.mat'),'allsubjs_Gradient','-v7.3')
disp('Gradient of all subjects [115,1] & group average have been saved!');
%result have been updated 1 feb


%% 
% LOAD SAVED DATA
% allsubjs_Gradient = load(fullfile(path_output,'Gradient_allsubjs.mat')).allsubjs_Gradient; % data of all HCP subjects
% Avg_Gradient = load(fullfile(path_output,'Gradient_GroupAvg.mat')).Avg_Gradient;

% %-------------------------------------------------------
%   use SVD to transfer each subject into avg group map
% %-------------------------------------------------------
allsubjs_Gradient_svd = nan(max(parc),length(subjs));

for s = 1:length(subjs)

    subj = num2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    temp = Avg_Gradient.*allsubjs_Gradient{s,1};
    temp_size = size(temp);

    %svd input must not contain nan
    % Get indices of finite values (non NaNs)
    [xf_l,~]=find(isfinite(temp));
    ix_finite = unique(xf_l);

    % Get indices of infinite values (NaNs)
    xi_l = setdiff([1:size(temp,1)]',xf_l);

    %exclude nans -->  vector without nans
    temp(unique(xi_l),:)=[];

    S = nan(temp_size);
    [S(ix_finite),~,T] = svd (temp,'econ');
    allsubjs_Gradient_svd(:,s) = allsubjs_Gradient{s,1}.*(T*S')';

end

%save data in output path
save(fullfile(path_output,'Gradient_allsubjs_to_avg.mat'),'allsubjs_Gradient_svd','-v7.3')
disp('Gradient of all subjects transfered to group average');


% Visualization in volume
% %---------------------------------------------------------------------------
% % we have a gradinet value for each brain region based on our parcellation (115)
% % set the gradient value to all voxels located in each specific region
% 
% parc_vol = parc_file.parc;
% U1_brain_vol = cell(length(subjs),1);
% 
% for s=1:length(subjs)
%     
%     region_val = U1_brain_svd{s,1};
%     vol =zeros(91,109,91);
%     
%     for i=1:max(parc1)       
%         idx = find((parc_vol ==i));
%         vol(idx) = region_val(i); 
%     end
%         U1_brain_vol{s,1} = vol;
% end
% 
% for i = 1:91; imagesc(vol(:, :, i)); axis equal; axis off;colorbar; caxis([-2, 2]); pause; end
% %---------------------------------------------------------------------------
% % plot corr(lh,rh) across all subjects
% Gradient_Hoacer=allsubjs_Gradient_svd;
% Gradient_Hoacer=G_HCP;
% %lh_label indx in parc label file --> 1:48,   97:103, 114
% %rh_label indx in parc label file --> 49:96, 104:110, 113
% [~,y_n]=find(isnan(Gradient_Hoacer)); %7 subjects has nan in some regions
% Gradient_Hoacer(:,unique(y_n))=[];
% 
% 
% lh_Gradient = [Gradient_Hoacer(1:48,:); Gradient_Hoacer(97:103,:); Gradient_Hoacer(114,:)];
% rh_Gradient = [Gradient_Hoacer(49:96,:);Gradient_Hoacer(104:110,:); Gradient_Hoacer(113,:)];
% 
% [x,y] = find(isnan(rh_Gradient));
% y=unique(y);
% 
% lh_Gradient(:,y)=[];
% rh_Gradient(:,y)=[];
% 
% 
% figure, plot(lh_Gradient, rh_Gradient, '.'); axis square, hold on
% xlabel("left hem."), ylabel("right hem."),
% title("Gradient, corr : "+corr(lh_Gradient(:), rh_Gradient(:)))





