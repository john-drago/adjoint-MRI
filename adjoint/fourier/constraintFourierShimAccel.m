function [ c, shimc ] = constraintFourierShimAccel( pSc, opt )

shimAccelConstr = opt.shimAccel_constr;

% Assume pSc is scaled to be between -1 and 1
pSc = pSc( : );

numCoils = opt.numZCoils;
numFourier = opt.numFourier_shim;
FBshimMat = opt.D2_shim;

shim_fourier_idx = opt.shim_idx;
shim_fourier_idx_rshp = reshape( shim_fourier_idx, [ numFourier, numCoils ] );
shim_fourier_sc = opt.scVec( shim_fourier_idx );
shim_fourier_sc_rshp = reshape( shim_fourier_sc, [ numFourier, numCoils ] );
shim_fourier = shim_fourier_sc .* pSc( shim_fourier_idx );
shim_fourier_rshp = reshape( shim_fourier, [ numFourier, numCoils ] );

shim = FBshimMat * shim_fourier_rshp;
shim_abs = abs( shim );

[ shim_max, shim_max_pos ] = max( shim_abs, [], 1 );

c_unsc = transpose( shim_max ) - shimAccelConstr;
c = c_unsc / shimAccelConstr;

if nargout > 1

    shim_max_vals_idx = sub2ind( size( shim_abs ),...
        shim_max_pos, 1:numCoils );

    shimc_unsc = zeros( opt.numVars, numCoils );

    dshimabsdshim = sign( shim( shim_max_vals_idx ) );
    dshimabsdshimf = transpose( dshimabsdshim ) .* FBshimMat( shim_max_pos, : );
    
    shimc_shimabs_shimf = sub2ind( size( shimc_unsc ),...
        transpose( shim_fourier_idx_rshp ),...
        repmat( transpose(1:numCoils), [ 1, numFourier ] ) );
    shimc_unsc( shimc_shimabs_shimf ) =...
        ( dshimabsdshimf ) .* transpose( shim_fourier_sc_rshp );

    shimc = shimc_unsc / shimAccelConstr;

end

end