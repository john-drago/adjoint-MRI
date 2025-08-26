function varArray = checkShimVariables( varArray, opt )
varNames = string( varArray( :, 1 ) );
if opt.numZCoils == 0
    shimVars = contains(varNames, "shim", 'ignorecase', true);
    varArray = varArray( ~shimVars, : );
end
end