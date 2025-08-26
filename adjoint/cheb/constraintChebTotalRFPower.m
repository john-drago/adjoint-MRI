function [ c, gradc ] = constraintChebTotalRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFtotalpowerconstr = opt.totalRFPower_constr;

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
c_breal_pwr = sum( breal_pwr, 2 );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numCheb, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numCheb, numXYCoils ] );
bimag_cheb = bimag_sc .* pSc( bimag_idx );
bimag_cheb_rshp = reshape( bimag_cheb, [ numCheb, numXYCoils ] );
bimag = TnT * bimag_cheb_rshp;
bimag_pwr = bimag.^2;
c_bimag_pwr = sum( bimag_pwr, 2 );

c_unsc = dutyCycle/( 2 * opt.Z0 * numTimePoints ) *...
    ( sum( c_breal_pwr + c_bimag_pwr ) ) - RFtotalpowerconstr;
c = c_unsc / RFtotalpowerconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Assign gradients with respect to breal (or bx)
    grad_breal = 2 * transpose( TnT ) * ( TnT * breal_cheb_rshp );
    gradc_unsc( breal_idx_rshp ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_breal ) .* breal_sc_rshp;
    
    % Assign gradients with respect to bimag (or by)
    grad_bimag = 2 * transpose( TnT ) * ( TnT * bimag_cheb_rshp );
    gradc_unsc( bimag_idx_rshp ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( grad_bimag ) .* bimag_sc_rshp;

    gradc = gradc_unsc / RFtotalpowerconstr;
end
end