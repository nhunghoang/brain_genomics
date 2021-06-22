%% calculating ALFF/fALFF
%%-- Neda SP, May 2021

%   Start with a timeseries before filtering T.
%	Then filter the signal in the range of 0.01-0,1Hz. Call this signal T_filtered.
%	ALFF is defined like this: rms(T_filtered)
%	fALFF is defined like this: rms(T_filtered)./rms(T)
%	rms is root mean square (~ standard deviation).

%------------data directory------------------
%platypus directory
path = fullfile('/data1/rubinov_lab/brain_genomics/data_HCP/');
timeseries_path = fullfile('/data1/rubinov_lab/brain_genomics/data_HCP/timeseries/');

%%my localcomputer directory
path = fullfile('/Users/neda/OneDrive - Vanderbilt/PHD/summer2021');
timeseries_path = fullfile('/Users/neda/OneDrive - Vanderbilt/PHD/summer2021/timeseries/');

%121 region names
region_label = string(readcell(fullfile(path, 'regions-121.txt')));


%get a list of 1206 subjects

Files=dir(fullfile(timeseries_path,'*.mat'));
Name = {Files.name};

%----------------------------ALFF/fALFF-----------------------------------%
TR = 0.720;        %Asampleperiod
LowCutoff = 0.01;
HighCutoff = 0.08;
sampleLength=1200;  %nDimTimePoints
sampleFreq 	 = 1/TR;

%load data : mat files
for sbj=1:size(Name,2) %loop through subjects
    
    temp_path = fullfile(timeseries_path, Name(1,sbj));
    p=Name(1,sbj);
    cd timeseries
    data = load(p{1});
    
    for ii=1:4        % avg between 4 fMRI scan matrices
        mat(ii,:,:) = cell2mat(timeseries(ii,1));
        timeseries = squeeze(mean(mat,1));  % we have tomeseries for 122 brain regions
    end
    
    
    timeseries = 2*abs(fft(timeseries))/sampleLength;  %122*1200
    
    
    %Get the frequency index
   
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
    
    
    %Calculating ALFF/fALFF
    ALFF = mean(timeseries(:,idx_LowCutoff:idx_HighCutoff),2);
    fALFF = sum(timeseries(:,idx_LowCutoff:idx_HighCutoff),2) ./ sum(timeseries(:,2:(paddedLength/2 + 1)),2);
    signal_properties(:,:,sbj) = [ALFF, fALFF];  % save ALFF/fALFF of all subjects: [121,2,#subjects]
    
end


