function pulse = checkShimConstraints(pulse, opt)
    % Get list of constraint names
    constraintList = keys(pulse.constraints);
    
    % Identify shim-related constraints
    shimVars = contains(constraintList, "shim", 'IgnoreCase', true);
    
    if opt.numZCoils == 0 && any(shimVars)
        % Remove shim-related constraints if no z coils available
        pulse.constraints(constraintList(shimVars)) = [];
    end
end