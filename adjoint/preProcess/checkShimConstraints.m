function pulse = checkShimConstraints( pulse, opt )

constraintList = keys( pulse.constraints );

if opt.numZCoils == 0
    shimVars = contains( constraintList, "shim", 'ignorecase', true);
    pulse.constraints( constraintList( shimVars ) ) = [];
end
end