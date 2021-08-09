function preprocess_hcp_volume(parc, path_output, clean_flag)
% preprocess_hcp(parc, path_output, path_input)
%   function inputs:
%     parc:         parcellation volume
%     path_output:  path to store output
%
%   function summary:
%     1. Load preprocessed volumetric RS-fMRI scans from HCP datasets
%     2. Parcellate brain voxel timeseries using input parcellation
%     2. Regress out movement, WM and CSF confounds, and filter

if ~exist('clean_flag', 'var')
    clean_flag = 1;
end

% set default paths
path_input = '/data1/datasets/hcp'; %default input directory
path_metadata = '/data1/rubinov_lab/datasets/hcp/hcp_prepro/fetch/metadata/';
path_subj = '/data1/rubinov_lab/brain_genomics/analyses_HCP/';
mkdir(path_output);

% get list of subjects and metadata
subj_list = h5read(fullfile(path_subj, 'subj_samp_assoc_order.hdf5'), '/subjects');
hl = load(fullfile(path_metadata, 'unrestricted_mrubinov_4_24_2019_15_25_31'));
metadata = struct('Subject', hl.data.Subject, 'QC_Issue', hl.data.QC_Issue);
metadata.Subject = string(metadata.Subject);

% get scans of interest
scan_name = "rfMRI_REST" + ["1_LR"; "1_RL"; "2_LR"; "2_RL"];
scan_tr = 0.72 * ones(size(scan_name));
scans = table(scan_name, scan_tr, 'VariableNames', {'Name', 'TR'});

% reshape and get parcells
parc = parc(:);
pmax = max(parc, [], 'all');

% loop over hcp subjects
parfor i = 1:length(subj_list)
    
    % get subject name and index
    subj = string(subj_list(i));
    disp(subj)
    ix_subj = find(metadata.Subject == subj);
    
    % make directory for each subj and open file for writing
    filename = fullfile(path_output, subj+".mat");
    
    if isfile(filename)
        continue
    else
        handle = matfile(filename, 'writable', true);
        handle.scans = scans;
        handle.qc = metadata.QC_Issue(ix_subj);
    end
    
    regressors = cell(height(scans), 1);
    Vp_clean = cell(height(scans), 1);
    for h = 1:height(scans)
        % check for quality control
        if ~isundefined(handle.qc)
            handle.errors(h, 1) = {'qc issue'};
            continue
        end
        
        scan_name = scans.Name{h};
        path_scan = fullfile(path_input, subj, 'MNINonLinear', 'Results', scan_name);
        
        % get regression data names
        rnams = [
            "Movement_Regressors", "Movement";
            scan_name+"_CSF", "CSF";
            scan_name+"_WM", "WM";
            "Movement_RelativeRMS", "Movt_RRMS"];
        
        % load and store regression data
        rdata = struct();
        for r = 1:size(rnams, 1)
            rdata.(rnams(r, 2)) = tryload(@load, fullfile(path_scan, rnams(r, 1)+".txt"));
        end
        regressors{h} = rdata;
        
        % make sure all regressors exist
        if any(structfun(@isempty, rdata))
            handle.errors(h, 1) = {'empty regressors'};
            continue
        end
        
        % get timeseries and brain mask
        V = tryload(@niftiread, fullfile(path_scan, scan_name+"_hp2000_clean.nii.gz"));
        
        % make sure data exists
        if isempty(V)
            handle.errors(h, 1) = {'empty data'};
            continue
        end
        tmax = size(V, 4);
        if tmax < 1200
            handle.errors(h, 1) = {'truncated data'};
            continue
        end
        V = double(reshape(V, [], tmax));
        
        %#ok<*AGROW> make regressors for model
        regr = rdata.Movement(:, 1:6);                      % motion
        regr = [regr [zeros(1, 6); regr(1:end-1,:)]];       % shifted
        regr = [regr rdata.CSF rdata.WM];                   % CSF, WM
        
        % make filter
        nyquist = (1 / scans.TR(h)) / 2;
        filt_kernel = fir1(ceil(tmax/6)*2-1, [0.01 nyquist-eps]/nyquist);
        
        % loop over parcels and apply regressors
        Vp_clean{h} = nan(pmax, tmax);
        for u = 1:pmax
            vp1 = mean(V(parc==u, :), 1, 'omitnan');
            if clean_flag
                Vp_clean{h}(u, :) = clean(vp1, regr, filt_kernel);
            else
                Vp_clean{h}(u, :) = vp1;
            end
        end
    end
    handle.Vp_clean = Vp_clean;
    handle.regressors = regressors;
end

end

function v = clean(v, regr, filt_kernel)

if ~isempty(v) && all(isfinite(v)) && isequal(size(regr, 1), numel(v))
    v = v(:);                               % reshape to column form
    m = mean(v);                            % get mean
    regr(:, end+1) = 1;                     % add intercept term
    v = v - regr * (regr \ v);              % remove regressors
    v = filtfilt(filt_kernel, 1, v);        % filter
    v = v - mean(v) + m;                    % restore mean
end

end

function data = tryload(load_function, full_filename)

[~, filename, ext] = fileparts(full_filename);

data = [];
if isfile(full_filename)
    for i = 1:10
        try
            temp_full_filename = tempname+"_"+filename+ext;     % generate random filename
            copyfile(full_filename, temp_full_filename);        % make a temporary copy of file
            data = load_function(temp_full_filename);           % load from a temporary copy
            delete(temp_full_filename)                          % delete temporary copy
            break
        catch
            disp("failed to load "+temp_full_filename+" on attempt "+i);
        end
    end
end

end
