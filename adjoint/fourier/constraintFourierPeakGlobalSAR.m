function [ c, gradc ] = constraintFourierPeakGlobalSAR( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numFourier = opt.numFourier_RF;
FBRFMat = opt.FBRF;
QGlobal = opt.QGlobal;

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
bcomp = transpose( FBRFMat * bcomp_rshp );
% breal = real( transpose( bcomp ) );
% bimag = imag( transpose( bcomp ) );

globalSARarr = permute( squeeze( real(...
    sum( conj(bcomp) .* ... 
    tensorprod( QGlobal, bcomp, 2, 1 ) , 1 ) ) ), [ 2 1 ] );

[ peakGlobalSAR, maxTimeidx ] = max( globalSARarr, [], 'all' );

c_unsc = ( peakGlobalSAR ) - peakGlobalSARconstr;
c = c_unsc / peakGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    gradf_zbar = 2 * QGlobal * ( bcomp( :, maxTimeidx ) );
    
    % Assign real and imaginary parts of the gradient
    grad_breal_brealf = real( FBRFMat( maxTimeidx, :) );
    grad_breal_bimagf = -imag( FBRFMat( maxTimeidx, :) );
    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +...
        transpose( real( gradf_zbar ) * ( grad_breal_brealf ) ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        transpose( real( gradf_zbar ) * ( grad_breal_bimagf ) ) .* bimag_fourier_sc_rshp;

    grad_bimag_brealf = imag( FBRFMat( maxTimeidx, :) );
    grad_bimag_bimagf = real( FBRFMat( maxTimeidx, :) );
    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +...
        transpose( imag( gradf_zbar ) * ( grad_bimag_brealf ) ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        transpose( imag( gradf_zbar ) * ( grad_bimag_bimagf ) ) .* bimag_fourier_sc_rshp;

    gradc = gradc_unsc / peakGlobalSARconstr;
end
end