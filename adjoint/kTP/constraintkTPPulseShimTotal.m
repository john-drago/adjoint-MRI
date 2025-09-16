function [ c, gradc ] = constraintkTPPulseShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

num_kTP = opt.num_kTP;
numZCoils = opt.numZCoils;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
shim_abs = abs( shim );

shim_sc_rshp = reshape( shim_sc, [ numZCoils, (num_kTP-1) ] );
shim_idx_rshp = reshape( shim_idx, [ numZCoils, (num_kTP-1) ] );
shim_rshp = reshape( shim, [ numZCoils, (num_kTP-1) ] );
shim_abs_rshp = reshape( shim_abs, [ numZCoils, (num_kTP-1) ] );

c_unsc = zeros( (num_kTP-1), 1 );

for kk = 1:(num_kTP-1)
    c_unsc( kk ) = -shimTotalConstr + sum( shim_abs_rshp( :, kk ) );
end

c = c_unsc / shimTotalConstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, (num_kTP-1) );
    
    for kk = 1:(num_kTP-1)
        gradc_unsc( shim_idx_rshp( :, kk ), kk ) = sign( shim_rshp( :, kk ) ) .* shim_sc_rshp( :, kk );
    end

    gradc = gradc_unsc / shimTotalConstr;

end

end