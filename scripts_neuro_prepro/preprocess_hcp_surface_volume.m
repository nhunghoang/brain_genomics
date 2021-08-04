function fx = preprocess_hcp_mmpmel(Parc, Parc_label)

%>>>>>>>>>>**********************************************<<<<<<<<<<<<<<<<
%   function inputs:
%                    Run >> parcellation_MMPmel.m 
%>>>>>>>>>>**********************************************<<<<<<<<<<<<<<<<
%** In summary, this function do this jobs:
% 1. Load preprocessed RS-fMRI scans from HCP datasets, both volumetric (.nii.gz) & CIFTI grayordinates data (.dtseries.nii)
% 2. Performs some extra preprocessing steps: filtering & regressing out, like:Movement, WM, CSF, offset
% 3. Parcellates brain voxel timeseries using Glasser(2016)_HCP_MMP1.0_5 + Melburne parcellations


% get paths to input and output
path_input = '/data1/datasets/hcp';                                                %read directory
path_accre = '/data1/rubinov_lab/brain_genomics/data_HCP/mmpmel/timeseries'; %write directory
mkdir(path_accre);

% get subject metadata
hl = load('/data1/rubinov_lab/datasets/hcp/hcp_prepro/fetch/metadata/unrestricted_mrubinov_4_24_2019_15_25_31');
metadata = struct('Subject', hl.data.Subject, 'QC_Issue', hl.data.QC_Issue);
metadata.Subject = string(metadata.Subject);

% get list of subject IDs to preprocess 
path_subj = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5';
subj_list = h5read(path_subj, '/subjects');

% get scans of interest
scans = table([], [], 'variablenames', {'Name', 'TR'});
scans = [scans;
    {'rfMRI_REST1_LR', 0.72}; {'rfMRI_REST1_RL', 0.72};
    {'rfMRI_REST2_LR', 0.72}; {'rfMRI_REST2_RL', 0.72}];

% parcellation indices, number of voxels, number of parcells, labels
Ix_parc = cellfun(@(x) find(x > 0), Parc, 'uniformoutput', false);
N_voxels = cellfun(@nnz, Ix_parc);
Pmax = cellfun(@(x) max(x(:)), Parc);
label = Parc_label;
%%
% offset = 60; % set an offset for data (remove first 60 timepoints in voxel timeseries) and regressors 

% loop over hcp subjects
parfor i = 1:length(subj_list)
    
    % get subject name and index
    subj = string(subj_list(i));
    ix = find(metadata.Subject == subj);
    assert(numel(ix) == 1, 'subject indexing error');
    
    % make directory for each subj and open file for writing
    mkdir(fullfile(path_accre, subj));
    filename = fullfile(path_accre, subj, 'processed_mmpmel.mat');
    
    % % temp speed-up %
    % vars = whos('-file', filename);
    % nams = cellfun(@string, {vars.name});
    % if any(nams == "V_clean") && any(nams == "regressors")
    %     fprintf('%s: already processed, skipping.\n', subj);
    %     continue;
    % end
    % % end temp speed-up %
    
    delete(filename)
    handle = matfile(filename, 'writable', true);
    handle.scans = scans;
    handle.qc = metadata.QC_Issue(ix);
    
    % check for quality control
    assert(isundefined(handle.qc), 'qc problem');
    % str = sprintf('%s: qc problem.', subj);
    % handle.err = str; disp(str)
    % continue
    % end
    
    regressors = cell(height(scans), 1);
    V_clean = cell(height(scans), 2);
    Vp_clean = cell(height(scans), 2);
    Vp_alff = cell(height(scans), 2);
    Vp_falff = cell(height(scans), 2);
    
    %loop through scans
    for h = 1:height(scans)
        scan_name = scans.Name{h};
        path2scan = fullfile(path_input, subj, 'MNINonLinear', 'Results', scan_name);
        
        % get regression data names
        rnams = {
            'Movement_Regressors', 'Movement';
            [scan_name '_CSF'], 'CSF';
            [scan_name '_WM'], 'WM';
            'Movement_RelativeRMS', 'Movt_RRMS'};
        
        % load and store regression data
        rdata = struct();
        for r = 1:size(rnams, 1)
            rdata.(rnams{r, 2}) = ...
                load(fullfile(path2scan, [rnams{r, 1} '.txt']));
        end
        regressors{h} = rdata;
        
        % make sure all regressors exist
        if any(structfun(@isempty, rdata))
            str = sprintf('%s: empty regressors.', subj);
            handle.(sprintf('scan_err%d', h)) = str; disp(str)
            continue
        end
        
        % exclude subjects with high motion
        if max(rdata.Movt_RRMS) > 0.2
            str = sprintf('%s_%d: motion problem: %f %f', subj, h, ...
                mean(rdata.Movt_RRMS), max(rdata.Movt_RRMS));
            handle.(sprintf('scan_err%d', h)) = str;
            disp(str)
        else
            fprintf('%s_%d: motion ok: %f %f\n', subj, h, ...
                mean(rdata.Movt_RRMS), max(rdata.Movt_RRMS));
        end
        
        %#ok<*AGROW> make regressors for model
        regr = rdata.Movement(:, 1:6);                      % motion
        regr = [regr [zeros(1, 6); regr(1:end-1,:)]];       % shifted
        regr = [regr rdata.CSF rdata.WM];                   % CSF, WM
        % regr = [power((1:size(regr,1)).', 0:2) regr];     % trends
        regr = regr(offset+1:end, :);                       % offset
        
        %loop through parcellations 
        % (Volume-based subcortical parcels & Surface-based Cortical parcels)
        for g = 1:numel(Parc)
            parc = Parc{g};
            ix_parc = Ix_parc{g};   % parcels' index from input
            n_voxels = N_voxels(g); % # of voxels in parcellation from input
            pmax = Pmax(g);         % # of parcels from input
            
            % get timeseries and brain mask
            switch g
                case 1        % for subcortical
                    file_name = [scan_name '_hp2000_clean.nii.gz'];
                    V = niftiread(fullfile(path2scan, file_name));
                    V = double(reshape(V, [], size(V, 4)));
                    
                case 2        % for cortical
                    file_name = [scan_name '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
                    cift = ft_read_cifti(fullfile(path2scan, file_name));
                    V = [cift.dtseries(cift.brainstructure==1, :);
                        cift.dtseries(cift.brainstructure==2, :)];
            end
            
            % V = V(:, offset+1:end);  %remove first 6o timepoints
            tmax = size(V, 2);
            
            % make filter
            nyquist = (1 / scans.TR(h)) / 2;
            % try
            filt_kernel = fir1(ceil(tmax/6)*2-1, [0.01 nyquist-eps]/nyquist);
            % catch
            % error(num2str([i h g]))
            % end
            
            Vp_clean{h, g} = zeros(pmax, tmax); 
            Vp_alff{h, g} = zeros(pmax, 1);
            Vp_falff{h, g} = zeros(pmax, 1);
            %------------------ parcellation --------------------------
            for t = 1:tmax
                vt = double(V(:, t)); %get voxel values at each timepoint
                Vp_clean{h, g}(:, t) = accumarray(parc(ix_parc), vt(ix_parc), [], @mean); %parcellated timeseries
            end
            
            %loop through parcellated timeseries
            for p = 1:pmax
                vp0 = Vp_clean{h, g}(p, :).';
                [vo, v] = clean(vp0, regr, filt_kernel);
                
                %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                %calculating ALFF/fALFF before filtering step
                nfft = 2^nextpow2(tmax); % compute power spectrum
                ft = fft(v, nfft) / tmax;
                freq = nyquist * linspace(0, 1, nfft/2+1).';
                powr = 2 * abs(ft(1:nfft/2+1, :));
                ALFF = sum(powr(freq > 0.01 & freq < 0.1));
                fALFF = Vp_alff{h, g}(p) / sum(powr);
                %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                
                Vp_clean{h, g}(p, :) = vo;
                Vp_alff{h, g}(p) = ALFF;
                Vp_falff{h, g}(p) = fALFF;
            end
            
            % get signal for voxels (voxels from parcels that we are interested in)
            
            V_clean{h, g} = zeros(n_voxels, tmax, 'single');            
            for j = 1:n_voxels  %loop through voxels
                v0 = double(V(ix_parc(j), :).');
                [vo, ~] = clean(v0, regr, filt_kernel); 
                V_clean{h, g}(j, :) = vo; 
            end
            
            
        end
    end
    
    
    handle.V_clean = V_clean;   %clean data for voxels (voxels of our parcellation)
    handle.Vp_clean = Vp_clean; %parcellated clean data
    handle.Vp_alff = Vp_alff;
    handle.Vp_falff = Vp_falff;
    handle.regressors = regressors;
end
end


            

function [vo,v] = clean(v, regr, filt_kernel)
% 1. Regressing out the confounds (Movement, CSF, WM)--> out : v signal
% 2. Filtering --> out : vo signal

if all(isfinite(v))
    %v = v(:);
    m = mean(v);
    v = v - regr * (regr \ v); % regressing out the confounds same as glmfit(regr, v) --> use v signal for ALFF calculation (non-filted data)
    vo = v;
    vo = filtfilt(filt_kernel, 1, vo);
    vo = vo - mean(vo) + m;
end
 
end
