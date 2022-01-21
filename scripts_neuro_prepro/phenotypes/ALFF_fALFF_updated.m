%% calculating ALFF/fALFF 
% Neda SP, Aug 2021 updated

% this script load timeseries matrices (.mat file) as input in forms of: [parcels, timepoints] or [allvoxels,timepoints]

% >>>>> processing steps:
%   >Start with a timeseries before filtering T.
%	>Then filter the signal in the range of 0.01 - 0,1 Hz. Call this signal T_filtered.
%	>ALFF is defined like this: rms(T_filtered)
%	>fALFF is defined like this: rms(T_filtered)./rms(T)
%	>rms is root mean square (~ standard deviation).


% get file info of valid subject scans 
ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5';%891
subjs = h5read(ts_order, '/subjects');

%------------timeseries directory------------------
%%--hoacer_sn_hth - new atlas with 115 parcels (filt & no-filt)

%hoacer_path = fullfile('/data1/rubinov_lab/brain_genomics/data_HCP_hoacer');
filt_hoacer_sn_hth_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth';
%nofilt_hoacer_sn_hth_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_nofilt';
hoacer_sn_hth_voxels = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_voxels';
mmpmel_filt = '/data1/rubinov_lab/brain_genomics/data_HCP/preprocessed_mmpmel/mmpmel_filt';


timeseries_path = fullfile(filt_hoacer_sn_hth_path, 'timeseries');
%get a list of timeseries
Files=dir(fullfile(timeseries_path,'*.mat'));
Name = {Files.name}';


%------- load parcellation ----------
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;

%loading mask for nonzero voxels
mask = load(fullfile(hoacer_sn_hth_voxels, 'mask')).mask;
mask = mask(:); 
mask_ind = find(mask);

%loop through subjects
%parcel
subj_ALFF = nan(max(parc(:)), length(subjs));
subj_fALFF = nan(max(parc(:)), length(subjs));
%voxel
subj_ALFF = nan(length(mask_ind(:)), length(subjs));
subj_fALFF = nan(length(mask_ind(:)), length(subjs));
next=1;

for s = 1:length(subjs)
    subj = int2str(subjs(s));
    disp([num2str(s), ': ', subj]);
    
    
    %load data : mat files
    mat_name = [subj '.mat'];
    %mat_name = ['processed_hoacer.mat'];
    temp_path = fullfile(timeseries_path, mat_name);
    data = load(temp_path);
    
    %data.Vp_clean --> scan metrices
    
    for ii=1:4        % loop through fMRI scan metrices
        
        timeseries = cell2mat(data.Vp_clean(ii,1));
        t = size(timeseries, 2);
        sampleLength=t; 
        
        if  isempty(timeseries)==1 %SKIP empty scans
            X = [num2str(ii),'th scan of',subj,'is empty!'];disp(X)
            timeseries =nan(max(parc(:)),sampleLength); %set the timeseries of empty scan-->nan
            No_scan(next,:) = cellstr(X);
            next=next+1;
            continue
        else    
            timeseries = double(2*abs(fft(timeseries))/sampleLength);  %powerspectrum
        end 
        
        %-----------------------------------------
        % ALFF for all parcels
        TR=0.7200; 
        %TR = data.scans(ii,2);%table
        LowCutoff = 0.01;
        HighCutoff = 0.08;  
        sampleFreq 	 = 1/TR;

        %-----Get the frequency index --------
        
        paddedLength = 2^nextpow2(sampleLength);
        if (LowCutoff >= sampleFreq/2) % All high included
            idx_LowCutoff = paddedLength/2 + 1;
        else % high cut off, such as freq > 0.01 Hz
            idx_LowCutoff = ceil(LowCutoff * paddedLength * TR + 1);
            % Change from round to ceil: idx_LowCutoff = round(LowCutoff *paddedLength *ASamplePeriod + 1);
        end
        if (HighCutoff>=sampleFreq/2)||(HighCutoff==0) % All low pass
            idx_HighCutoff = paddedLength/2 + 1;
        else % Low pass, such as freq < 0.08 Hz
            idx_HighCutoff = fix(HighCutoff *paddedLength *TR + 1);
            % Change from round to fix: idx_HighCutoff	=round(HighCutoff *paddedLength *ASamplePeriod + 1);
        end
        %-----------------------------------
        
        %Calculating ALFF/fALFF for each fMRI scan
        disp([num2str(ii), 'th scan of ', subj,' :']);
        ALFF(:,ii) = mean(timeseries(:,idx_LowCutoff:idx_HighCutoff),2);disp(['ALFF completed']);
        fALFF(:,ii) = sum(timeseries(:,idx_LowCutoff:idx_HighCutoff),2) ./ sum(timeseries(:,2:(paddedLength/2 + 1)),2);disp(['fALFF completed']);
   
    end % through scans
    
    % average the subject's ALFF/fALFF across their scans
    avg_subj_ALFF = nanmean(ALFF, 2);
    avg_subj_fALFF = nanmean(fALFF, 2);%mean of non nans
    
    subj_ALFF(:,s) = avg_subj_ALFF;
    subj_fALFF(:,s) = avg_subj_fALFF;
    disp([subj, ' : scan averaging finished']);
  

end %through subjects


%% save result matrix
cd /data1/rubinov_lab/Neda/ALFF
%parcel-based alff
save('ALFF_hoacer_sn_hth_nofilt.mat','subj_ALFF','-v7.3') % ALFF matrix 
save('fALFF_hoacer_sn_hth_nofilt.mat','subj_fALFF','-v7.3') % fALFF matrix 
%voxel-based ALFF
save('ALFF_hoacer_sn_hth_voxels.mat','subj_ALFF','-v7.3') % ALFF matrix 
save('fALFF_hoacer_sn_hth_voxels.mat','subj_fALFF','-v7.3') % fALFF matrix 

%% loading results

%parcel data
subj_ALFF_parcel = load('ALFF_hoacer_sn_hth_filt.mat').subj_ALFF;
subj_fALFF_parcel = load('fALFF_hoacer_sn_hth_filt.mat').subj_fALFF;

%voxel data
subj_ALFF_voxels = load('ALFF_hoacer_sn_hth_voxels.mat').subj_ALFF;
subj_fALFF_voxels = load('fALFF_hoacer_sn_hth_voxels.mat').subj_fALFF;

%% -------parcellation------------------
%----------------------------------------
   parc_mask=parc(mask_ind);     
   for s=1:length(subjs)+1     
         avg_subj_ALFF=subj_ALFF_voxels(:,s);
         avg_subj_fALFF=subj_fALFF_voxels(:,s);
        
        for j = 1:max(parc_mask(:))
            ALFF_115parc(j,s) = mean(avg_subj_ALFF(parc_mask==j),1);
            fALFF_115parc(j,s) = mean(avg_subj_fALFF(parc_mask==j),1);
        end
   end


%parceled alff
save('ALFF_hoacer_sn_hth_vox_parcelled.mat','ALFF_115parc','-v7.3') % ALFF matrix 
save('fALFF_hoacer_sn_hth_vox_parcelled.mat','fALFF_115parc','-v7.3') % fALFF matrix 

%% -------------------------ploting------------------------------------------

%115 parcel
%lh_label indx in parc label file --> 1:48,   97:103, 114
%rh_label indx in parc label file --> 49:96, 104:110, 113

figure,plot(ALFF_115parc(1:48,1),subj_ALFF_parcel(1:48,1),'.');axis square,hold on,
xlabel("parcellated voxel timeseries"),ylabel("parcel timeseries"),
title("hoacer-sn-hth, ALFF for parcels 1:48")

%------------------correlation across subjects
%ALFF
subj_ALFF_filt = ALFF_115parc;

lh_ALFF_filt = [subj_ALFF_filt(1:48,:);subj_ALFF_filt(97:103,:); subj_ALFF_filt(114,:)];
rh_ALFF_filt = [subj_ALFF_filt(49:96,:);subj_ALFF_filt(104:110,:); subj_ALFF_filt(113,:)];

figure, plot(lh_ALFF_filt, rh_ALFF_filt, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("hoacer-sn-hth:  filt-ALFF corr(lh,rh):"+corr(lh_ALFF_filt(:), rh_ALFF_filt(:)))

%fALFF
lh_fALFF = [subj_fALFF(1:50,:);subj_fALFF(109:116,:); subj_fALFF(120,:)];
rh_fALFF = [subj_fALFF(51:100,:);subj_fALFF(101:108,:); subj_fALFF(121,:)];
for h=1:size(lh_fALFF,1)

    f_crr_mat(h,:)=corr(lh_fALFF(h,:)', rh_fALFF(h,:)');
end

figure, plot(lh_fALFF, rh_fALFF, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("re-preprocessed data:  fALFF")



%% -----------------within one subject across ROIs
%ALFF
lh_ALFF_filt = [subj_ALFF(2:50,1);subj_ALFF(109:116,1); subj_ALFF(120,1)];
rh_ALFF_filt = [subj_ALFF(52:100,1);subj_ALFF(101:108,1); subj_ALFF(121,1)];

figure, plot(lh_ALFF_filt, rh_ALFF_filt, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("re-preproc data, one subj,  ALFF corr : "+corr(lh_ALFF_filt, rh_ALFF_filt))

%fALFF
lh_fALFF = [subj_fALFF(2:50,1);subj_fALFF(109:116,1); subj_fALFF(120,1)];
rh_fALFF = [subj_fALFF(52:100,1);subj_fALFF(101:108,1); subj_fALFF(121,1)];

figure, plot(lh_fALFF, rh_fALFF, '.'); axis square, hold on
xlabel("left hem."), ylabel("right hem."),
title("re-preproc data, one subj,  fALFF corr : "+corr(lh_fALFF, rh_fALFF))