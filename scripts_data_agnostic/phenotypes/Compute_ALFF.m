

function [ALFF_parcelvalues, ALFF_voxelvalues, fALFF_parcelvalues, fALFF_voxelvalues] = Compute_ALFF(project, path_output)

%-----------------------------------------------------------------
%
%   function summary:
%     this script load timeseries matrices (.mat file) as input in forms of: [parcels, timepoints] or [allvoxels,timepoints]
%       >>>>> processing steps:
%        >Start with a timeseries before filtering T.
%	     >Then filter the signal in the range of 0.01 - 0,1 Hz. Call this signal T_filtered.
%	     >ALFF is defined like this: rms(T_filtered)
%	     >fALFF is defined like this: rms(T_filtered)./rms(T)
%	     >rms is root mean square (~ standard deviation).
%
%   function inputs:
%     project:      label: human connectome project("hcp") or UK-biobank("ukb")
%     path_output:  path to store output Like: path_output = '/data1/rubinov_lab/Neda/ALFF'
%
%   function outouts:
%     ALFF_parcelvalues, fALFF_parcelvalues :  contain ALFF and fALFF result matrices for 115 brain regions based on "hoacer_sn_hth" parcellation file :[115, number of subjects]
%     ALFF_voxelvalues,  fALFF_voxelvalue   :  contain ALFF and fALFF result matrices for all brain voxels: [number of voxels, number of subjects]
%
%
%  The result matrices for "HCP" and "UKB" are saved in the following paths:
%---------UKB------------- 
%  path_output = '/data1/rubinov_lab/Neda/ALFF_UKB/'
%   >> ALFF_all_voxels.mat; 
%   >> fALFF_all_voxels.mat;
%   >> ALFF_hoacer_sn_hth_vox_parcelled.mat;
%   >> fALFF_hoacer_sn_hth_vox_parcelled.mat

%---------HCP-----------
%  path_output = '/data1/rubinov_lab/Neda/ALFF/'
%   >> ALFF_all_voxels.mat; 
%   >> fALFF_all_voxels.mat;
%   >> ALFF_hoacer_sn_hth_vox_parcelled.mat;
%   >> fALFF_hoacer_sn_hth_vox_parcelled.mat
%
%*Note: the calculation pipeline for both cohorts is the same!
%
% Neda Sardaripour, 2021 
%-----------------------------------------------------------------
%-----------------------------------------------------------------

%----- load parcellation: hoacer_sn_hth (115 brain regions) ------
parc_file = load('/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/parcellation.mat');
parc = parc_file.parc;
parc = parc(:);
parc_label = parc_file.name;
        
%loading mask for nonzero voxels
% mask = load(fullfile(hoacer_sn_hth_voxels, 'mask')).mask;
mask = nonzeros(parc);
mask_ind = find(mask);


%frequency band for calculating ALFF
 LowCutoff = 0.01;
 HighCutoff = 0.08;

%------------------------------------------------------------------------------------------       
%>>>>> switch between two projects: 1. Human Connectome Project(hcp) 2.UKBiobank(ukb) <<<<<%
%------------------------------------------------------------------------------------------
switch project
    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>Human Connectome Project (HCP)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------
    case "hcp"

        disp(['ALFF/fALFF calculation in HCP dataset has started...']);
        % input path
        hoacer_sn_hth_voxels = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth_voxels/timeseries';
        %Old_hoacer_sn_hth_voxels = '/data1/rubinov_lab/brain_genomics/data_HCP/2021_hoacer_sn_hth_voxels';
        %path output: (from function input)
        if ~exist(path_output)
             mkdir(path_output);
        end
        
        % get file info of valid subject scans 
        ts_order = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5';%890
        subjs = h5read(ts_order, '/subjects');  

%         %get a list of timeseries
%         Files=dir(fullfile(hoacer_sn_hth_voxels,'*.mat'));
%         Name = {Files.name}';
%         Name = Name(1:890,1);

        %output variables for voxelbased calculation
        subj_ALFF = nan(length(mask_ind(:)), length(subjs));
        subj_fALFF = nan(length(mask_ind(:)), length(subjs));
        next=1; % for storing no data subjects

        
        %loop through subjects
        for s = 1:length(subjs)
            subj = int2str(subjs(s));
            disp([num2str(s), ': ', subj]);
            
            
            %load data : mat files
            mat_name = [subj '.mat'];
            %mat_name = ['processed_hoacer.mat'];
            temp_path = fullfile(hoacer_sn_hth_voxels, mat_name);
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
            avg_subj_ALFF = mean(ALFF, 2, 'omitnan');
            avg_subj_fALFF = mean(fALFF, 2, 'omitnan');
            
            subj_ALFF(:,s) = avg_subj_ALFF;
            subj_fALFF(:,s) = avg_subj_fALFF;
            disp([subj, ' : calculation finished']);
          
        
        end %through subjects

        %>>>>>>>output variable of function<<<<<<<
        ALFF_voxelvalues = subj_ALFF;
        fALFF_voxelvalues = subj_fALFF;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/ALFF_UKB

        %voxel-based ALFF
        save(fullfile(path_output,'ALFF_all_voxels.mat'),'subj_ALFF','No_scan','-v7.3') % ALFF matrix 
        save(fullfile(path_output,'fALFF_all_voxels.mat'),'subj_fALFF','No_scan','-v7.3') % fALFF matrix 
        disp(['voxel based result has saved in output_path']);

    %--------------------------------------------------------------------------------------
    %       >>>>>>>>>>>>>>>>>>>>UK Biobank (UKB)<<<<<<<<<<<<<<<<<<<<<
    %--------------------------------------------------------------------------------------    
    %Feb 15, 2021 updated With new final_subjs.txt order file

    case "ukb"
        
        disp(['ALFF/fALFF calculation in ukb dataset has started...']);
        % input path
        path_input = '/data1/rubinov_lab/brain_genomics/data_UKB/hoacer_sn_hth_voxels/timeseries';
        timeseries_path = path_input;
        %path_output = '/data1/rubinov_lab/Neda/ALFF_UKB/'
        
        if ~exist(path_output)
             mkdir(path_output);
        end
                
        % get list of subjects
        UKB_subjs_order = fullfile('/data1/rubinov_lab/Neda/UKB_subjs_list/');
        subjs = load(fullfile(UKB_subjs_order, 'final_subjs.txt'));
 

        %get scan of interest
        scan_name = "rfMRI";
        scan_tr = 0.735 * ones(size(scan_name));
        scans = table(scan_name, scan_tr, 'VariableNames', {'Name', 'TR'});
        
        %loop through subjects
        
        for s = 1:length(subjs)
            subj = num2str(subjs(s));
            disp([num2str(s), ': ', subj]);
            
            %load data : mat files
            mat_name = [subj];
            
            temp_path = fullfile(timeseries_path, mat_name);
            data = load(temp_path);
            
            %data.Vp_clean --> scan metrices
                            
            %output variables for voxelbased calculation
            subj_ALFF = nan(length(mask_ind(:)), length(subjs));
            subj_fALFF = nan(length(mask_ind(:)), length(subjs));
            next=1; % for storing no data subjects 

        %    for ii=1:4       % loop through fMRI scan metrices: JUST ONE SCAN FOR UKB
                ii=1;
                timeseries = cell2mat(data.Vp_clean(ii,1));
                t = size(timeseries, 2);
                sampleLength=t; 
                
                if  isempty(timeseries)==1 %SKIP empty scans
                    X = [num2str(ii),'th scan of',subj,'is empty!'];disp(X)
                    timeseries =nan(length(mask_ind(:)),490); %set the timeseries of empty scan-->nan
                    No_scan(next,:) = cellstr(X);
                    next=next+1;
                else    
                    timeseries = double(2*abs(fft(timeseries))/sampleLength);  %powerspectrum
                end 
                
                %-----------------------------------------
                % ALFF for all parcels
                TR=data.scans.tr;  
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
           
         %   end %through scans
            
            % average the subject's ALFF/fALFF across their scans
            avg_subj_ALFF = mean(ALFF, 2, 'omitnan');
            avg_subj_fALFF = mean(fALFF, 2, 'omitnan');
            
            subj_ALFF(:,s) = avg_subj_ALFF;
            subj_fALFF(:,s) = avg_subj_fALFF;
            disp([subj, ' : calculation finished']);
          
        
        end %through subjects
        
        %>>>>>>>output variable of function<<<<<<<
        ALFF_voxelvalues = subj_ALFF;
        fALFF_voxelvalues = subj_fALFF;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        %-------save result matrix in "path_output" ----------
        %cd /data1/rubinov_lab/Neda/ALFF_UKB

        %voxel-based ALFF
        save(fullfile(path_output,'ALFF_all_voxels.mat'),'subj_ALFF','No_scan','-v7.3') % ALFF matrix 
        save(fullfile(path_output,'fALFF_all_voxels.mat'),'subj_fALFF','No_scan','-v7.3') % fALFF matrix 
        disp(['voxel based result has saved in output path']);

end



%  >> After calculating and saving voxel based ALFF/fALFF in HCP or UKB cohorts:
%  >> calculate values for each brain parcels (hoacer_sn_hth parcellation file)
%
% %% loading results
% %voxel data
% subj_ALFF_voxels = load('ALFF_hoacer_sn_hth_voxels.mat').subj_ALFF;
% subj_fALFF_voxels = load('fALFF_hoacer_sn_hth_voxels.mat').subj_fALFF;

%--------------------------------------------------------------------------------
% -------parcellation: caluculate parcel based ALFF/fALFF values------------------
%--------------------------------------------------------------------------------
disp(['parcellation has started...']);

% get mask indices
idx_parc = find(parc);
nmax = nnz(parc);
mask = parc(idx_parc);

  % parcellation
  for s=1:length(subjs)        
        for j = 1:max(mask(:))          
            ALFF_115parc(j,s)= nanmean(nonzeros(subj_ALFF(mask==j,s)));
            fALFF_115parc(j,s)= nanmean(nonzeros(subj_fALFF(mask==j,s)));
        end
  end 

        
 %>>>>>>>>output variable of function<<<<<<<<<<<<<<
  ALFF_parcelvalues = ALFF_115parc;
  fALFF_parcelvalues =  fALFF_115parc;
 %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% no scan
if exist(No_scan)
      for k=1:length(No_scan)
          temp = char(No_scan(k,1));
          no_data(k,:) = temp(12:22);
      end 
else
      no_data=[];
    
end

%----------save data in path_output----------
%parcel based alff
save(fullfile(path_output, 'ALFF_hoacer_sn_hth_vox_parcelled.mat'),'ALFF_115parc','no_data','-v7.3') 
save(fullfile(path_output,'fALFF_hoacer_sn_hth_vox_parcelled.mat'),'fALFF_115parc','no_data','-v7.3')  
disp(['parcellated results has saved in output path']);

end




