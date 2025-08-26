function [ c, gradc ] = constraintMPpTxPulseMaxRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFmaxpowerconstr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
num_ORSP = opt.num_ORSP;
num_MPSP = opt.num_MPSP;

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );
brealORSP_rshp = reshape( brealORSP, [ num_ORSP, numXYCoils ] );
brealORSP_pwr_rshp = brealORSP_rshp.^2;
c_brealORSP_pwr = transpose( sum( brealORSP_pwr_rshp, 1 ) );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );
brealMPSP_rshp = reshape( brealMPSP, [ num_MPSP, numXYCoils ] );
brealMPSP_pwr_rshp = brealMPSP_rshp.^2;
c_brealMPSP_pwr = transpose( sum( brealMPSP_pwr_rshp, 1 ) );

c_breal_pwr = c_brealORSP_pwr + c_brealMPSP_pwr;

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );
bimagORSP_rshp = reshape( bimagORSP, [ num_ORSP, numXYCoils ] );
bimagORSP_pwr_rshp = bimagORSP_rshp.^2;
c_bimagORSP_pwr = transpose( sum( bimagORSP_pwr_rshp, 1 ) );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );
bimagMPSP_rshp = reshape( bimagMPSP, [ num_MPSP, numXYCoils ] );
bimagMPSP_pwr_rshp = bimagMPSP_rshp.^2;
c_bimagMPSP_pwr = transpose( sum( bimagMPSP_pwr_rshp, 1 ) );

c_bimag_pwr = c_bimagORSP_pwr + c_bimagMPSP_pwr;

channelPowersUncorrected = c_breal_pwr + c_bimag_pwr;
[ ~, highPowerChannel ] = max( channelPowersUncorrected );

c_unsc = dutyCycle/( 2 * opt.Z0 * numTimePoints ) *...
    ( channelPowersUncorrected( highPowerChannel ) ) - RFmaxpowerconstr;
c = c_unsc / RFmaxpowerconstr;

if nargout > 1
    
    % Initialize gradient vector
    % gradc = sparse( opt.numVars, 1 );
    gradc_unsc = zeros( opt.numVars, 1 );

    % high power channel indices
    highPowerChannelORSPIdx = ( ( ( highPowerChannel - 1) * num_ORSP ) + ( 1:num_ORSP ) ).';
    highPowerChannelMPSPIdx = ( ( ( highPowerChannel - 1) * num_MPSP ) + ( 1:num_MPSP ) ).';

    % Assign gradients with respect to breal (or bx)
    gradc_unsc( brealORSP_idx( highPowerChannelORSPIdx ), 1 ) = ( ( 2 * dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        brealORSP( highPowerChannelORSPIdx ) .* brealORSP_sc( highPowerChannelORSPIdx );
    gradc_unsc( brealMPSP_idx( highPowerChannelMPSPIdx ), 1 ) = ( ( 2 * dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        brealMPSP( highPowerChannelMPSPIdx ) .* brealMPSP_sc( highPowerChannelMPSPIdx );
    
    % Assign gradients with respect to bimag (or by)
    gradc_unsc( bimagORSP_idx( highPowerChannelORSPIdx ), 1 ) = ( ( 2 * dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        bimagORSP( highPowerChannelORSPIdx ) .* bimagORSP_sc( highPowerChannelORSPIdx );
    gradc_unsc( bimagMPSP_idx( highPowerChannelMPSPIdx ), 1 ) = ( ( 2 * dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        bimagMPSP( highPowerChannelMPSPIdx ) .* bimagMPSP_sc( highPowerChannelMPSPIdx );

    gradc = gradc_unsc / RFmaxpowerconstr;
end
end