function varArray = checkShimVariables(varArray, opt)
    % Extract variable names as string array
    varNames = string(varArray(:,1));
    
    % Identify shim-related variables
    shimVars = contains(varNames, "shim", 'IgnoreCase', true);
    
    if opt.numZCoils == 0
        % Remove all shim-related variables if no z coils available
        varArray = varArray(~shimVars, :);
        
    elseif any(shimVars)
        % If there are shim variables, check their associated values
        shimValues = cell2mat(varArray(shimVars, 2));
        
        % If all shim-related values are zero, remove them
        if all(shimValues == 0)
            varArray = varArray(~shimVars, :);
        end
    end
end
