function [ c, gradc ] = constraintSpokesMaxRFPower( pSc, opt, spokes )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFmaxpowerconstr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
pulseLength = spokes.pulseLength;
numSpokesCurr = length( opt.breal_idx ) / opt.numXYCoils;
numXYCoils = opt.numXYCoils;

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numXYCoils, numSpokesCurr ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numXYCoils, numSpokesCurr ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, numSpokesCurr ] );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, numSpokesCurr ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, numSpokesCurr ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, numSpokesCurr ] );

c_bimag_pwr = zeros( numXYCoils, 1 );
c_breal_pwr = zeros( numXYCoils, 1 );
RF_spoke_pwr = zeros( numSpokesCurr, 1 );

for nn = 1:numSpokesCurr
    
    spokeIdx = spokes.numSpokes - nn + 1;
    colIdx = numSpokesCurr - nn + 1;

    RF_spoke_int_pwr = spokes.dt * sum( spokes.RF_spokes{ spokeIdx }.^2 );
    RF_spoke_pwr( colIdx ) = RF_spoke_int_pwr;
    
    c_breal_pwr_spoke = ( breal_rshp( :, colIdx ).^2 * RF_spoke_int_pwr );
    c_bimag_pwr_spoke = ( bimag_rshp( :, colIdx ).^2 * RF_spoke_int_pwr );

    c_breal_pwr = c_breal_pwr + c_breal_pwr_spoke;
    c_bimag_pwr = c_bimag_pwr + c_bimag_pwr_spoke;
end

channelPowersUncorrected = c_breal_pwr + c_bimag_pwr;
[ highChannelPower, highPowerChannel ] = max( channelPowersUncorrected );

maxpwr = dutyCycle / ( 2 * opt.Z0 * pulseLength ) *...
    ( highChannelPower  );

c_unsc = maxpwr - RFmaxpowerconstr;
c = c_unsc / RFmaxpowerconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % high power channel indices
    % highPowerChannelIdx = ( ( ( highPowerChannel - 1) * numTimePoints ) + ( 1:numTimePoints ) ).';

    % Assign gradients with respect to breal (or bx)
    grad_breal = 2 * ( transpose( breal_rshp( highPowerChannel, : ) ) .* RF_spoke_pwr );
    gradc_unsc( breal_idx_rshp( highPowerChannel, : ), 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * pulseLength ) ) *...
        ( grad_breal ) .* transpose( breal_sc_rshp( highPowerChannel, : ) );

    % Assign gradients with respect to bimag (or by)
    grad_bimag = 2 * ( transpose( bimag_rshp( highPowerChannel, : ) ) .* RF_spoke_pwr );
    gradc_unsc( bimag_idx_rshp( highPowerChannel, : ), 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * pulseLength ) ) *...
        ( grad_bimag ) .* transpose( bimag_sc_rshp( highPowerChannel, : ) );

    gradc = gradc_unsc / RFmaxpowerconstr;
end

end