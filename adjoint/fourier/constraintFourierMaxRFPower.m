function [ c, gradc ] = constraintFourierMaxRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFmaxpowerconstr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
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
breal_pwr = breal.^2;
c_breal_pwr = transpose( sum( breal_pwr, 1 ) );

bimag = imag( bcomp );
bimag_pwr = bimag.^2;
c_bimag_pwr = transpose( sum( bimag_pwr, 1 ) );

channelPowersUncorrected = c_breal_pwr + c_bimag_pwr;
[ ~, highPowerChannel ] = max( channelPowersUncorrected );

maxpwr = dutyCycle / ( 2 * opt.Z0 * numTimePoints ) *...
    ( channelPowersUncorrected( highPowerChannel )  );

c_unsc = maxpwr - RFmaxpowerconstr;
c = c_unsc / RFmaxpowerconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Assign gradients with respect to breal (or bx)
    grad_breal_brealf = 2 * transpose( real( FBRFMat ) ) * ( breal( :, highPowerChannel ) );
    grad_breal_bimagf = 2 * transpose( -imag( FBRFMat ) ) * ( breal( :, highPowerChannel ) );

    gradc_unsc( breal_fourier_idx_rshp( :, highPowerChannel ) ) =...
        gradc_unsc( breal_fourier_idx_rshp( :, highPowerChannel ) ) +...
        ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal_brealf ) .* breal_fourier_sc_rshp( :, highPowerChannel );

    gradc_unsc( bimag_fourier_idx_rshp( :, highPowerChannel ) ) =...
        gradc_unsc( bimag_fourier_idx_rshp( :, highPowerChannel ) ) +...
        ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal_bimagf ) .* bimag_fourier_sc_rshp( :, highPowerChannel );
    
    % Assign gradients with respect to bimag (or by)
    grad_bimag_brealf = 2 * transpose( imag( FBRFMat ) ) * ( bimag( :, highPowerChannel ) );
    grad_bimag_bimagf = 2 * transpose( real( FBRFMat ) ) * ( bimag( :, highPowerChannel ) );

    gradc_unsc( breal_fourier_idx_rshp( :, highPowerChannel ) ) =...
        gradc_unsc( breal_fourier_idx_rshp( :, highPowerChannel ) ) +...
        ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag_brealf ) .* breal_fourier_sc_rshp( :, highPowerChannel );

    gradc_unsc( bimag_fourier_idx_rshp( :, highPowerChannel ) ) =...
        gradc_unsc( bimag_fourier_idx_rshp( :, highPowerChannel ) ) +...
        ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag_bimagf ) .* bimag_fourier_sc_rshp( :, highPowerChannel );

    gradc = gradc_unsc / RFmaxpowerconstr;
end
end