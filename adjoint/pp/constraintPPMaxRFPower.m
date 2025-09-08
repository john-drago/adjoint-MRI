function [ c, gradc ] = constraintPPMaxRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFmaxpowerconstr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
varsToTimepoints_RF = opt.varsToTimepoints_RF;

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numVarsPerChannel_RF, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numVarsPerChannel_RF, numXYCoils ] );
breal_shape = breal_sc .* pSc( breal_idx );
breal_shape_rshp = reshape( breal_shape, [ numVarsPerChannel_RF, numXYCoils ] );
breal = varsToTimepoints_RF * breal_shape_rshp;
breal_pwr = breal.^2;
c_breal_pwr = transpose( sum( breal_pwr, 1 ) );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_shape = bimag_sc .* pSc( bimag_idx );
bimag_shape_rshp = reshape( bimag_shape, [ numVarsPerChannel_RF, numXYCoils ] );
bimag = varsToTimepoints_RF * bimag_shape_rshp;
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
    grad_breal = 2 * transpose( varsToTimepoints_RF ) * ( varsToTimepoints_RF * breal_shape_rshp( :, highPowerChannel ) );
    gradc_unsc( breal_idx_rshp( :, highPowerChannel ) ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal ) .* breal_sc_rshp( :, highPowerChannel );

    % Assign gradients with respect to bimag (or by)
    grad_bimag = 2 * transpose( varsToTimepoints_RF ) * ( varsToTimepoints_RF * bimag_shape_rshp( :, highPowerChannel ) );
    gradc_unsc( bimag_idx_rshp( :, highPowerChannel ) ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag ) .* bimag_sc_rshp( :, highPowerChannel );

    gradc = gradc_unsc / RFmaxpowerconstr;
end

end