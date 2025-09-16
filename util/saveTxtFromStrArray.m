function saveTxtFromStrArray(filename, strArray)
% This function saves the strings in a string array to the file
% specified by filename.
%
% Between every entry of the strArray, the function adds '\n' to force
% a new line in the .txt file.

% Open the file for writing in text mode
fid = fopen(filename, 'wt');

if fid == -1
    error('Failed to open file: %s', filename);
end

try
    % Iterate over each string in the array
    for nn = 1:numel(strArray)
        % Retrieve the current string
        currentStr = strArray(nn);
        
        % Write the string as-is without escaping
        fprintf(fid, '%s\n', currentStr);
    end
catch ME
    % Close the file in case of an error
    fclose(fid);
    rethrow(ME);
end

% Close the file after writing
fclose(fid);
end
