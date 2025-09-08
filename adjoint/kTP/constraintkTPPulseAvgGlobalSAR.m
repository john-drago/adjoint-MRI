function [ c, gradc ] = constraintkTPPulseAvgGlobalSAR( pSc, opt )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;
dtvec = transpose( opt.dtvec( : ) );
pulseLength = opt.pulseLength;

dt_RF_noslew = dtvec( opt.RF_idx );

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
% breal_idx_rshp = reshape( breal_idx, [ numXYCoils, num_kTP ] );
% breal_sc_rshp = reshape( breal_sc, [ numXYCoils, num_kTP ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, num_kTP ] );
breal_half_rshp = breal_rshp / 2;

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
% bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, num_kTP ] );
% bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, num_kTP ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
bimag_half_rshp = bimag_rshp / 2;

bcomp_rshp = transpose( complex( breal_rshp, bimag_rshp ) );
% bcomp_rshp_t = transpose( bcomp_rshp );
bcomp_half_rshp = transpose( complex( breal_half_rshp, bimag_half_rshp ) );

b_pwr_high = ( ctranspose( bcomp_rshp ) .* dt_RF_noslew ) * bcomp_rshp;
b_pwr_half = ( bcomp_half_rshp' .* (2 * opt.RFSlewTime ) ) * bcomp_half_rshp;

b_pwr = b_pwr_high + b_pwr_half;

avgGlobalSAR = real( squeeze( sum( QGlobal .* b_pwr, [ 1 2 ] ) ) );

c_unsc = ( (dutyCycle) / pulseLength ) * ( avgGlobalSAR ) - avgGlobalSARconstr;

c = c_unsc / avgGlobalSARconstr;
% c = c_unsc;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    gradf_zbar_b_rshp = ( ( (dutyCycle) / pulseLength ) * ( 2 ) )...
        * ( ( bcomp_rshp .* ( transpose( dt_RF_noslew ) + (0.5 * opt.RFSlewTime ) ) ) * conj( QGlobal ) );
    gradf_zbar_b = reshape( transpose( gradf_zbar_b_rshp ), [ num_kTP * numXYCoils ,1 ] );

    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx ) = real( gradf_zbar_b ) .* breal_sc;
    gradc_unsc( bimag_idx ) = imag( gradf_zbar_b ) .* bimag_sc;

    % assign to gradc
    gradc = gradc_unsc / avgGlobalSARconstr;
    % gradc = gradc_unsc;

end
end