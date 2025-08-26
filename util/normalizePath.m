function normalizedPath = normalizePath(inputPath)
    % normalizePath: Resolves '.', '..', and '~' in a file path.
    % Inputs:
    %   inputPath - The path to normalize (can be relative or absolute)
    % Outputs:
    %   normalizedPath - The resolved absolute path without '.', '..', '~'

    % Expand '~' to home directory on both Unix and Windows
    if startsWith(inputPath, '~')
        if ispc
            homeDir = getenv('USERPROFILE'); % Windows home directory
        else
            homeDir = getenv('HOME'); % Unix home directory
        end
        if strcmp(inputPath, '~')
            inputPath = homeDir;
        else
            inputPath = fullfile(homeDir, inputPath(2:end));
        end
    end

    % Convert the path to an absolute path if it's not already
    if ~isabs(inputPath)
        inputPath = fullfile(pwd, inputPath);
    end

    % Split the path into parts
    pathParts = strsplit(inputPath, filesep);

    % Initialize a stack for the normalized path
    stack = {};

    for ii = 1:length(pathParts)
        part = pathParts{ii};
        if strcmp(part, '.') || isempty(part)
            % Skip current directory references
            continue;
        elseif strcmp(part, '..')
            % Remove the last directory from the stack if possible
            if ~isempty(stack)
                stack(end) = [];
            else
                error('Invalid path: Too many ".." references.');
            end
        else
            % Add the directory to the stack
            stack{end + 1} = part; %#ok<AGROW>
        end
    end

    % Join the stack into the normalized path
    normalizedPath = strjoin(stack, filesep);

    % Prepend the root if it's an absolute path
    if ispc && startsWith(inputPath, '\') % Windows root-relative path
        normalizedPath = ['\' normalizedPath];
    elseif ispc && length(inputPath) >= 2 && inputPath(2) == ':' % Windows drive letter
        normalizedPath = [pathParts{1} filesep normalizedPath];
    elseif startsWith(inputPath, filesep) % Unix absolute path
        normalizedPath = [filesep normalizedPath];
    end
end

%% Helper Functions
% ----------------------------------------------------------------------- %
function tf = isabs(path)
    % isabs: Checks if a path is absolute.
    % Works for both Windows and Unix-like paths.
    if ispc
        tf = ~isempty(regexp(path, '^[A-Za-z]:\\', 'once')); % Windows absolute path
    else
        tf = startsWith(path, '/'); % Unix absolute path
    end
end
% ----------------------------------------------------------------------- %