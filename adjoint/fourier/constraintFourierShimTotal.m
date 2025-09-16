function [ c, gradc ] = constraintFourierShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSc is scaled to be between -1 and 1
pSc = pSc( : );


numCoils = opt.numZCoils;
numFourier = opt.numFourier_shim;
FBshimMat = opt.FBshim;

shim_fourier_idx = opt.shim_idx;
shim_fourier_idx_rshp = reshape( shim_fourier_idx, [ numFourier, numCoils ] );
shim_fourier_sc = opt.scVec( shim_fourier_idx );
shim_fourier_sc_rshp = reshape( shim_fourier_sc, [ numFourier, numCoils ] );
shim_fourier = shim_fourier_sc .* pSc( shim_fourier_idx );
shim_fourier_rshp = reshape( shim_fourier, [ numFourier, numCoils ] );

shim = FBshimMat * shim_fourier_rshp;
shim_abs = abs( shim );
shim_abs_sum = sum( shim_abs, 2 );
[ shim_abs_sum_max, shim_abs_sum_max_idx ] = max( shim_abs_sum );

c_unsc = shim_abs_sum_max - shimTotalConstr;

c = c_unsc / shimTotalConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 1 );

    dshimabsdshim = sign( shim( shim_abs_sum_max_idx, : ) );
    dshimabsdshimf = transpose( dshimabsdshim ) .* FBshimMat( shim_abs_sum_max_idx, : );

    gradc_unsc( transpose( shim_fourier_idx_rshp ) ) =...
        ( dshimabsdshimf ) .* transpose( shim_fourier_sc_rshp );

    gradc = gradc_unsc / shimTotalConstr;

end

end