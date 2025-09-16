function [ c, gradc ] = constraintFourierAvgGlobalSAR( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
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
bcomp = ( FBRFMat * bcomp_rshp );

b_pwr = ctranspose( bcomp ) * bcomp;

avgGlobalSAR = real( sum( QGlobal .* b_pwr, [ 1 2 ] ) );

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgGlobalSAR ) - avgGlobalSARconstr;
c = c_unsc / avgGlobalSARconstr;

if nargout > 1

    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )

    % Determine contribution to average global SAR from each point
    gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
        * ( bcomp * conj( squeeze( QGlobal ) ) );

    % Assign real and imaginary parts of the gradient
    grad_breal_brealf = real( FBRFMat );
    grad_breal_bimagf = -imag( FBRFMat );
    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +....
        ( transpose( grad_breal_brealf ) * real( gradf_zbar_rshp ) ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        ( transpose( grad_breal_bimagf ) * real( gradf_zbar_rshp ) ) .* bimag_fourier_sc_rshp;

    grad_bimag_brealf = imag( FBRFMat );
    grad_bimag_bimagf = real( FBRFMat );
    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +....
        ( transpose( grad_bimag_brealf ) * imag( gradf_zbar_rshp ) ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        ( transpose( grad_bimag_bimagf ) * imag( gradf_zbar_rshp ) ) .* bimag_fourier_sc_rshp;

    gradc = gradc_unsc / avgGlobalSARconstr;
    
end
end