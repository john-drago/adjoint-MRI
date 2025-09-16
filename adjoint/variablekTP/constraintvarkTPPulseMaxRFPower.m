function [ c, gradc ] = constraintvarkTPPulseMaxRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

maxRFPower_constr = opt.maxRFPower_constr;

dutyCycle = opt.dutyCycle;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;

dt_idx = opt.dt_idx;
dt_sc = opt.scVec( dt_idx );
pulseLength = sum( dt_sc .* pSc( dt_idx ) );

dt_blip_idx = dt_idx( 2:2:end );
dt_blip_sc = opt.scVec( dt_blip_idx );

dt_RF_idx = dt_idx( 1:2:end );
dt_RF_sc = opt.scVec( dt_RF_idx );
dt_RF = ( dt_RF_sc .* pSc( dt_RF_idx ) ).';
dt_RF_noslew = dt_RF - 2 * opt.RFSlewTime;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal_idx_rshp = reshape( breal_idx, [ numXYCoils, num_kTP ] );
breal_sc_rshp = reshape( breal_sc, [ numXYCoils, num_kTP ] );

breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, num_kTP ] );
breal_half_rshp = breal_rshp / 2;
breal_rshp_pwr = breal_rshp.^2;
breal_half_rshp_pwr = breal_half_rshp.^2;

breal_rshp_pwr_sc = sum( breal_rshp_pwr .* dt_RF_noslew, 2 );
breal_half_rshp_pwr_sc = sum( breal_half_rshp_pwr * (2 * opt.RFSlewTime ), 2);

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, num_kTP ] );
bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, num_kTP ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
bimag_half_rshp = bimag_rshp / 2;
bimag_rshp_pwr = bimag_rshp.^2;
bimag_half_rshp_pwr = bimag_half_rshp.^2;

bimag_rshp_pwr_sc = sum( bimag_rshp_pwr .* dt_RF_noslew, 2);
bimag_half_rshp_pwr_sc = sum( bimag_half_rshp_pwr * (2 * opt.RFSlewTime ), 2);

channelPowersUncorrected =...
    breal_rshp_pwr_sc + breal_half_rshp_pwr_sc +...
    bimag_rshp_pwr_sc + bimag_half_rshp_pwr_sc;

[ highPower, highPowerChannel ] = max( channelPowersUncorrected );

maxpwr = dutyCycle/( 2 * opt.Z0 ) * ( highPower  );

c_unsc = maxpwr - maxRFPower_constr * pulseLength;

% c = c_unsc / maxRFPower_constr;
c = c_unsc;

if nargout > 1
    
     % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Assign gradients with respect to breal (or bx)
    gradc_breal_unsc = (2 * opt.RFSlewTime) * ( breal_rshp( highPowerChannel, : ) / 2 ) + 2 * ( dt_RF_noslew .* breal_rshp( highPowerChannel, : ) );
    gradc_breal_sc = ( ( dutyCycle ) / ( 2 * opt.Z0 ) ) * gradc_breal_unsc .* breal_sc_rshp( highPowerChannel, : );
    gradc_unsc( breal_idx_rshp( highPowerChannel, : ), 1 ) = gradc_breal_sc;
    
    % Assign gradients with respect to bimag (or by)
    gradc_bimag_unsc = (2 * opt.RFSlewTime) * ( bimag_rshp( highPowerChannel, : ) / 2 ) + 2 * dt_RF_noslew .* bimag_rshp( highPowerChannel, : );
    gradc_bimag_sc = ( ( dutyCycle ) / ( 2 * opt.Z0 ) ) * gradc_bimag_unsc .* bimag_sc_rshp( highPowerChannel, : );
    gradc_unsc( bimag_idx_rshp( highPowerChannel, : ), 1 ) = gradc_bimag_sc;

    % Assign gradients with respect to dt
    gradc_dt_unsc = ( ( ( dutyCycle ) / ( 2 * opt.Z0 ) ) * ( breal_rshp_pwr( highPowerChannel, : ) + bimag_rshp_pwr( highPowerChannel, : ) )...
        .* dt_RF_sc.' ) -...
        maxRFPower_constr * dt_RF_sc.';
    gradc_unsc( dt_RF_idx, 1 ) = gradc_dt_unsc;

    gradc_unsc( dt_blip_idx, 1) = - maxRFPower_constr * dt_blip_sc;

    % gradc = gradc_unsc / maxRFPower_constr;
    gradc = gradc_unsc;

end

end