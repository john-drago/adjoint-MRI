function [ c, gradc ] = constraintPWCPulseTotalRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFtotalpowerconstr = opt.totalRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_pwr = breal.^2;
breal_pwr_rshp = reshape( breal_pwr, [ numTimePoints, numXYCoils ] );
c_breal_pwr = sum( breal_pwr_rshp, 2 );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_pwr = bimag.^2;
bimag_pwr_rshp = reshape( bimag_pwr, [ numTimePoints, numXYCoils ] );
c_bimag_pwr = sum( bimag_pwr_rshp, 2 );

c_unsc = dutyCycle/( 2 * opt.Z0 * numTimePoints ) *...
    ( sum( c_breal_pwr + c_bimag_pwr ) ) - RFtotalpowerconstr;
c = c_unsc / RFtotalpowerconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Assign gradients with respect to breal (or bx)
    grad_breal = 2 * breal;
    gradc_unsc( breal_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal ) .* breal_sc;
    
    % Assign gradients with respect to bimag (or by)
    grad_bimag = 2 * bimag;
    gradc_unsc( bimag_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag ) .* bimag_sc;

    gradc = gradc_unsc / RFtotalpowerconstr;
end
end