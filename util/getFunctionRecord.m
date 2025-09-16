function evalfnRecord = getFunctionRecord(currFile)
% getFunctionRecord - Reads the content of a file line by line and stores
% each line as an element in a string array.
%
% Inputs:
%   currFile - The name of the file to read (string or char).
%
% Outputs:
%   evalfnRecord - A string array containing each line of the file.

% Ensure the file exists
if ~isfile(currFile)
    error('File "%s" does not exist.', currFile);
end

% Open the file for reading
fid = fopen(currFile, 'r');
if fid == -1
    error('Failed to open file: %s', currFile);
end

% Initialize the record as an empty string array
evalfnRecord = string.empty;

try
    % Read the file line by line
    while true
        tline = fgetl(fid);
        if ~ischar(tline) % End of file
            break;
        end
        % Append the line to the record
        evalfnRecord((end+1), 1) = string(tline);
    end
catch ME
    % Ensure the file is closed in case of an error
    fclose(fid);
    rethrow(ME);
end

% Close the file after reading
fclose(fid);

end
