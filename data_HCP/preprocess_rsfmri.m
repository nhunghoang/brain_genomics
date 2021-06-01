function preprocessing_hcp(path_output, parc, typ)

% nhung path modification
path_output = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_seq';
log_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/prepro_logs_seq';

% number of parcels
parc = parc(:);
pmax = max(parc);

% get paths to input and output

% get subject metadata
path_code = fileparts(mfilename('fullpath'));
hdl = load([path_code '/../fetch/metadata/unrestricted_mrubinov_4_24_2019_15_25_31']);
metadata = struct('Subject', hdl.data.Subject, 'QC_Issue', hdl.data.QC_Issue);

% get scans of interest
scans = cell2table({ ...
    'rfMRI_REST1_LR', 0.72;
    'rfMRI_REST1_RL', 0.72;
    'rfMRI_REST2_LR', 0.72;
    'rfMRI_REST2_RL', 0.72}, ...
    'variablenames', {'Name', 'TR'});

% get subject directories
path_input = '/data1/datasets/hcp';
dsubj = dir(path_input);
dsubj = {dsubj.name};
dsubj = dsubj(cellfun(@(n) ~isempty(regexp(n, '[0-9]*')), dsubj)); %#ok<RGXP1>

%%

% nhung: only loop over subjects with missing ts 
% mfile = '/data1/rubinov_lab/brain_genomics/data_HCP/subjs_missing_ts.txt';
% fid = fopen(mfile);
% data = textscan(fid, '%s');
% fclose(fid);
% subj_list = string(data{:});

% loop over subjects
% for i = 1:length(subj_list)
for i = 1:length(dsubj)
    % get subject name and index
    subj_i = dsubj{i};
    %subj_i = subj_list{i};
    disp([num2str(i), ': ', subj_i]);
    
    % open file for writing
    warning off MATLAB:DELETE:FileNotFound
    filename = fullfile(path_output, [subj_i '.mat']);
    
    % nhung 
    %if isfile(filename)
    %    disp([num2str(i), ': ', subj_i, '(done)']);
    %    continue 
    %end 
    %disp([num2str(i), ': ', subj_i]);
    logID = fopen(fullfile(log_dir, [subj_i '.log']), 'w');
    
    delete(filename)
    hdl = matfile(filename);
    hdl.scans = scans;
    hdl.errors = cell(height(scans), 1);
    hdl.warnings = cell(height(scans), 1);
    hdl.qc = metadata.QC_Issue(metadata.Subject==str2num(subj_i));
    
    regressors = cell(height(scans), 1);
    Vp = cell(height(scans), 1);
    for h = 1:height(scans)
        % check for quality control
        if ~isundefined(hdl.qc)
            hdl.errors(h, 1) = {'qc issue'};
            fprintf(logID, '%s (%d): qc_issue\n', subj_i, h); % nhung
            continue
        end
        
        scan_h = scans.Name{h};
        path_scan = fullfile(path_input, subj_i, 'MNINonLinear', 'Results', scan_h);
        
        % get regression data names
        rnams = {
            'Movement_Regressors', 'Movement';
            [scan_h '_CSF'], 'CSF';
            [scan_h '_WM'], 'WM';
            'Movement_RelativeRMS', 'Movt_RRMS'};
        
        % load and store regression data
        try
            rdata = struct();
            for r = 1:size(rnams, 1)
                rdata.(rnams{r, 2}) = dlmread(fullfile(path_scan, [rnams{r, 1} '.txt']));
            end
        catch
            hdl.errors(h, 1) = {'no regressors'};
            fprintf(logID, '%s (%d): no regressors\n', subj_i, h); % nhung
            continue
        end
        regressors{h} = rdata;
        
        % check for high motion
        if max(rdata.Movt_RRMS) > 0.2
            hdl.warnings(h, 1) = {'high movement'};
        else
            hdl.warnings(h, 1) = {'no warnings'};
        end
        
        %#ok<*AGROW> make regressors for model
        %regr = rdata.Movement(:, 1:6);                      % motion
        %regr = [regr [zeros(1, 6); regr(1:end-1,:)]];       % shifted
        %regr = [regr rdata.CSF rdata.WM];                   % CSF, WM
        
        % get timeseries and brain mask
        switch typ
            case 'surface'
                name_scan = [scan_h '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
                try
                    cift = ft_read_cifti(fullfile(path_scan, name_scan));
                catch
                    hdl.errors(h, 1) = {'no data'};
                    continue
                end
                V = [cift.dtseries(cift.brainstructure==1, :);
                    cift.dtseries(cift.brainstructure==2, :)];
                t = size(V, 2);
            case 'volume'
                name_scan = [scan_h '_hp2000_clean.nii.gz'];
                for j = 1:5 % 5 attempts to get the data
                    try
                        V = niftiread(fullfile(path_scan, name_scan));
                        hdl.errors(h, 1) = {'no error'};
                        fprintf(logID, '%s (%d): no error\n', subj_i, h); % nhung
                        break;
                    catch ME
                        fprintf(logID, '%s (%d): %s\n', subj_i, h, ME.identifier); % nhung
                        hdl.errors(h, 1) = {'no data'};
                        pause(10*j)
                        continue
                    end
                end
                if strcmp(hdl.errors(h, 1), 'no data')
                    continue;
                end
                t = size(V, 4);
                V = double(reshape(V, [], t));
        end
        
        % make filter
        % nyquist = (1 / scans.TR(h)) / 2;
        % filt_kernel = fir1(ceil(t/6)*2-1, [0.01 nyquist-eps]/nyquist);
        
        % loop over parcellations
        if isempty(parc)
            Vp{h} = V;
        else
            Vp{h} = nan(pmax, t, 'single');
            for u = 1:pmax
                Vp{h}(u, :) = mean(V(parc==u, :), 1, 'double', 'omitnan');
            end
        end
        
    end
    
    hdl.timeseries = Vp;
    hdl.regressors = regressors;
    fclose(logID); % nhung 
end

end

function v = clean(v, regr, filt_kernel)

if all(isfinite(v))
    v = v(:);
    m = mean(v);
    v = v - regr * (regr \ v);
    v = filtfilt(filt_kernel, 1, v);
    v = v - mean(v) + m;
end

end
