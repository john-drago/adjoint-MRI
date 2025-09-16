function [ c, gradc ] = constraintChebMaxRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFmaxpowerconstr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numCheb = opt.numCheb_RF;
TnT = opt.Tn( :, 1:numCheb );

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numCheb, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numCheb, numXYCoils ] );
breal_cheb = breal_sc .* pSc( breal_idx );
breal_cheb_rshp = reshape( breal_cheb, [ numCheb, numXYCoils ] );
breal = TnT * breal_cheb_rshp;
breal_pwr = breal.^2;
c_breal_pwr = transpose( sum( breal_pwr, 1 ) );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numCheb, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numCheb, numXYCoils ] );
bimag_cheb = bimag_sc .* pSc( bimag_idx );
bimag_cheb_rshp = reshape( bimag_cheb, [ numCheb, numXYCoils ] );
bimag = TnT * bimag_cheb_rshp;
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
    grad_breal = 2 * TnT' * ( TnT * breal_cheb_rshp( :, highPowerChannel ) );
    gradc_unsc( breal_idx_rshp( :, highPowerChannel ) ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal ) .* breal_sc_rshp( :, highPowerChannel );

    % Assign gradients with respect to bimag (or by)
    grad_bimag = 2 * TnT' * ( TnT * bimag_cheb_rshp( :, highPowerChannel ) );
    gradc_unsc( bimag_idx_rshp( :, highPowerChannel ) ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag ) .* bimag_sc_rshp( :, highPowerChannel );

    gradc = gradc_unsc / RFmaxpowerconstr;
end

end