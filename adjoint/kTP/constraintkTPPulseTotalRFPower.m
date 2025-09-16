function [ c, gradc ] = constraintkTPPulseTotalRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

totalRFPower_constr = opt.totalRFPower_constr;

dutyCycle = opt.dutyCycle;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
dtvec = transpose( opt.dtvec( : ) );
pulseLength = opt.pulseLength;

dt_RF_noslew = dtvec( opt.RF_idx );

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
% breal_idx_rshp = reshape( breal_idx, [ numXYCoils, num_kTP ] );
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
% bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, num_kTP ] );
bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, num_kTP ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
bimag_half_rshp = bimag_rshp / 2;
bimag_rshp_pwr = bimag_rshp.^2;
bimag_half_rshp_pwr = bimag_half_rshp.^2;

bimag_rshp_pwr_sc = sum( bimag_rshp_pwr .* dt_RF_noslew, 2);
bimag_half_rshp_pwr_sc = sum( bimag_half_rshp_pwr * (2 * opt.RFSlewTime ), 2);


c_unsc = dutyCycle / ( 2 * opt.Z0 * pulseLength ) *...
    ( sum( breal_rshp_pwr_sc + breal_half_rshp_pwr_sc + bimag_rshp_pwr_sc + bimag_half_rshp_pwr_sc ) )...
    - totalRFPower_constr;

c = c_unsc / totalRFPower_constr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Assign gradients with respect to breal (or bx)
    gradc_breal_unsc = (2 * opt.RFSlewTime) * ( breal_rshp / 2 ) + 2 * ( dt_RF_noslew .* breal_rshp );
    gradc_breal_sc_rshp = ( ( dutyCycle ) / ( 2 * opt.Z0 * pulseLength ) ) * gradc_breal_unsc .* breal_sc_rshp;
    gradc_breal_sc = reshape( gradc_breal_sc_rshp, [ numXYCoils*num_kTP, 1 ] );
    gradc_unsc( breal_idx, 1 ) = gradc_breal_sc;
    
    % Assign gradients with respect to bimag (or by)
    gradc_bimag_unsc = (2 * opt.RFSlewTime) * (bimag_rshp / 2 ) + 2 * dt_RF_noslew .* bimag_rshp;
    gradc_bimag_sc_rshp = ( ( dutyCycle ) / ( 2 * opt.Z0 * pulseLength ) ) * gradc_bimag_unsc .* bimag_sc_rshp;
    gradc_bimag_sc = reshape( gradc_bimag_sc_rshp, [ numXYCoils*num_kTP, 1 ] );
    gradc_unsc( bimag_idx, 1 ) = gradc_bimag_sc;

    gradc = gradc_unsc / totalRFPower_constr;
    % gradc = gradc_unsc;
end
end