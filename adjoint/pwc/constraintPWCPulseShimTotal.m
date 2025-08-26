function [ c, gradc ] = constraintPWCPulseShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numTimePoints = opt.numTimePoints;
numZCoils = opt.numZCoils;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
shim_abs = abs( shim );

shim_sc_rshp = reshape( shim_sc, [ numTimePoints, numZCoils ] );
shim_idx_rshp = reshape( shim_idx, [ numTimePoints, numZCoils ] );
shim_rshp = reshape( shim, [ numTimePoints, numZCoils ] );
shim_abs_rshp = reshape( shim_abs, [ numTimePoints, numZCoils ] );

c_unsc = sum( shim_abs_rshp, 2 ) - shimTotalConstr;

c = c_unsc / shimTotalConstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numTimePoints );

    time_idx_rshp = repmat( transpose( 1:numTimePoints), [ 1 numZCoils ] );
    gradc_rshp = ( sign( shim_rshp ) .*  shim_sc_rshp );
    
    lin_idx = sub2ind( size( gradc_unsc ), shim_idx_rshp, time_idx_rshp );

    gradc_unsc( lin_idx ) = gradc_rshp;

    gradc = gradc_unsc / shimTotalConstr;

end

end