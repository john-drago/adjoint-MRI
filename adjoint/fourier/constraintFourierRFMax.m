function [ c, gradc ] = constraintFourierRFMax( pSc, opt )

RFMaxConstr = opt.RFMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numXYCoils = opt.numXYCoils;
numFourier = opt.numFourier_RF;
FBRFMat = opt.FBRF;

breal_fourier_idx = opt.breal_idx;
breal_fourier_idx_rshp = reshape( breal_fourier_idx, [ numFourier, numXYCoils ] );
breal_fourier_sc = opt.scVec( breal_fourier_idx );
breal_fourier_sc_rshp = reshape( breal_fourier_sc, [ numFourier, numXYCoils ] );
breal_fourier = breal_fourier_sc .* pSc( breal_fourier_idx );
breal_fourier_rshp = reshape( breal_fourier, [ numFourier, numXYCoils ] );

bimag_fourier_idx = opt.bimag_idx;
bimag_fourier_idx_rshp = reshape( bimag_fourier_idx, [ numFourier, numXYCoils ] );
bimag_fourier_sc = opt.scVec( bimag_fourier_idx );
bimag_fourier_sc_rshp = reshape( bimag_fourier_sc, [ numFourier, numXYCoils ] );
bimag_fourier = bimag_fourier_sc .* pSc( bimag_fourier_idx );
bimag_fourier_rshp = reshape( bimag_fourier, [ numFourier, numXYCoils ] );

bcomp_rshp = complex( breal_fourier_rshp, bimag_fourier_rshp );
bcomp = FBRFMat * bcomp_rshp;

breal = real( bcomp );
bimag = imag( bcomp );

breal_abs = abs( breal );
bimag_abs = abs( bimag );

[ breal_max, breal_max_pos ] = max( breal_abs, [], 1 );
[ bimag_max, bimag_max_pos ] = max( bimag_abs, [], 1 );

c_unsc = [ transpose( breal_max ); transpose( bimag_max ) ] - RFMaxConstr;
c = c_unsc / RFMaxConstr;

if nargout > 1

    breal_max_vals_idx = sub2ind( size( breal_abs ),...
        breal_max_pos, 1:numXYCoils );
    bimag_max_vals_idx = sub2ind( size( bimag_abs ),...
        bimag_max_pos, 1:numXYCoils );

    gradc_unsc = zeros( opt.numVars, 2*numXYCoils );
    
    dbrealdbrealf = +real( FBRFMat( breal_max_pos, : ) );
    dbrealdbimagf = -imag( FBRFMat( breal_max_pos, : ) );
    dbimagdbrealf = +imag( FBRFMat( bimag_max_pos, : ) );
    dbimagdbimagf = +real( FBRFMat( bimag_max_pos, : ) );

    dbrealabsdbreal = sign( breal( breal_max_vals_idx ) );
    dbimagabsdbimag = sign( bimag( bimag_max_vals_idx ) );
    
    dbrealabsdbrealf = transpose( dbrealabsdbreal ) .* dbrealdbrealf;
    dbrealabsdbimagf = transpose( dbrealabsdbreal ) .* dbrealdbimagf;

    dbimagabsdbrealf = transpose( dbimagabsdbimag ) .* dbimagdbrealf;
    dbimagabsdbimagf = transpose( dbimagabsdbimag ) .* dbimagdbimagf;
    
    gradc_brealabs_brealf = sub2ind( size( gradc_unsc ),...
        transpose( breal_fourier_idx_rshp ),...
        repmat( transpose(1:numXYCoils), [ 1, numFourier ] ) );

    gradc_brealabs_bimagf = sub2ind( size( gradc_unsc ),...
        transpose( bimag_fourier_idx_rshp ),...
        repmat( transpose(1:numXYCoils), [ 1, numFourier ] ) );

    gradc_unsc( gradc_brealabs_brealf ) =...
        ( dbrealabsdbrealf ) .* transpose( breal_fourier_sc_rshp );
    gradc_unsc( gradc_brealabs_bimagf ) =...
        ( dbrealabsdbimagf ) .* transpose( bimag_fourier_sc_rshp );
    
    gradc_bimagabs_brealf = sub2ind( size( gradc_unsc ),...
        transpose( breal_fourier_idx_rshp ),...
        repmat( numXYCoils + transpose(1:numXYCoils), [ 1, numFourier ] ) );

    gradc_bimagabs_bimagf = sub2ind( size( gradc_unsc ),...
        transpose( bimag_fourier_idx_rshp ),...
        repmat( numXYCoils + transpose(1:numXYCoils), [ 1, numFourier ] ) );

    gradc_unsc( gradc_bimagabs_brealf ) =...
        ( dbimagabsdbrealf ) .* transpose( breal_fourier_sc_rshp );
    gradc_unsc( gradc_bimagabs_bimagf ) =...
        ( dbimagabsdbimagf ) .* transpose( bimag_fourier_sc_rshp );

    gradc = gradc_unsc / RFMaxConstr;

end

end