function data = load_data(load_function, full_filename)

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
