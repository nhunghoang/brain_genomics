function v = clean_timeseries(v, regr, filt_kernel)

if ~isempty(v) && all(isfinite(v)) && isequal(size(regr, 1), numel(v))
    v = v(:);                               % reshape to column form
    m = mean(v);                            % get mean
    regr(:, end+1) = 1;                     % add intercept term
    v = v - regr * (regr \ v);              % remove regressors
    v = filtfilt(filt_kernel, 1, v);        % filter
    v = v - mean(v) + m;                    % restore mean
end

end
