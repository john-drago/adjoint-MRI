function [ c, gradc ] = constraintFourierPeakLocalSAR( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numFourier = opt.numFourier_RF;
FBRFMat = opt.FBRF;
VOPs = opt.VOPs;

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
bcomp_t = transpose( bcomp );
% breal = real( transpose( bcomp ) );
% bimag = imag( transpose( bcomp ) );

if ~ismatrix( VOPs )

    localSARarr = permute( squeeze( real(...
        sum( conj(bcomp) .* ...
        permute(...
        tensorprod( VOPs, bcomp, 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) ) ), [ 2 1 ] );

else

    localSARarr = transpose( real( ( VOPs * conj( bcomp_t ) ) .* bcomp_t ) );

end

[ peaklSAR, peaklSAR_li ] = max( localSARarr, [], 'all' );

c_unsc = ( peaklSAR ) - peakLocalSARconstr;
c = c_unsc / peakLocalSARconstr;

if nargout > 1

    % Determine peak VOP and time point
    [ highPowerVOP, maxTimeidx ] = ind2sub( size(localSARarr), peaklSAR_li );

    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    if ~ismatrix( VOPs )
        gradf_zbar = 2 * squeeze( VOPs( :, :, highPowerVOP ) ) * ( bcomp( :, maxTimeidx ) );
    else
        gradf_zbar = 2 * VOPs * ( bcomp( :, maxTimeidx ) );
    end

    % Assign real and imaginary parts of the gradient
    grad_breal_brealf = real( FBRFMat( maxTimeidx, :) );
    grad_breal_bimagf = -imag( FBRFMat( maxTimeidx, :) );

    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +...
        transpose( real( gradf_zbar ) * grad_breal_brealf ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        transpose( real( gradf_zbar ) * grad_breal_bimagf ) .* bimag_fourier_sc_rshp;

    grad_bimag_brealf = imag( FBRFMat( maxTimeidx, :) );
    grad_bimag_bimagf = real( FBRFMat( maxTimeidx, :) );

    gradc_unsc( breal_fourier_idx_rshp ) = gradc_unsc( breal_fourier_idx_rshp ) +...
        transpose( imag( gradf_zbar ) * grad_bimag_brealf ) .* breal_fourier_sc_rshp;
    gradc_unsc( bimag_fourier_idx_rshp ) = gradc_unsc( bimag_fourier_idx_rshp ) +...
        transpose( imag( gradf_zbar ) * grad_bimag_bimagf ) .* bimag_fourier_sc_rshp;

    gradc = gradc_unsc / peakLocalSARconstr;

end
end