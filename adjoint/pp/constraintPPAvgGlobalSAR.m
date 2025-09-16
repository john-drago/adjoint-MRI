function [ c, gradc ] = constraintPPAvgGlobalSAR( pSc, opt )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
varsToTimepoints_RF = opt.varsToTimepoints_RF;
QGlobal = opt.QGlobal;

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numVarsPerChannel_RF, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numVarsPerChannel_RF, numXYCoils ] );
breal_shape = breal_sc .* pSc( breal_idx );
breal_shape_rshp = reshape( breal_shape, [ numVarsPerChannel_RF, numXYCoils ] );
breal_rshp = varsToTimepoints_RF * breal_shape_rshp;

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_cheb = bimag_sc .* pSc( bimag_idx );
bimag_cheb_rshp = reshape( bimag_cheb, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_rshp = varsToTimepoints_RF * bimag_cheb_rshp;

bcomp_rshp = complex( breal_rshp, bimag_rshp );

b_pwr = ctranspose( bcomp_rshp ) * bcomp_rshp;

avgGlobalSAR = real( sum( QGlobal .* b_pwr, [ 1 2 ] ) );

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgGlobalSAR ) - avgGlobalSARconstr;
c = c_unsc / avgGlobalSARconstr;

if nargout > 1

    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )

    % Determine contribution to average global SAR from each point
    gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
        * ( bcomp_rshp * conj( squeeze( QGlobal ) ) );

    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx_rshp ) = ( transpose( varsToTimepoints_RF ) * real( gradf_zbar_rshp ) )...
        .* breal_sc_rshp;
    gradc_unsc( bimag_idx_rshp ) = ( transpose( varsToTimepoints_RF ) * imag( gradf_zbar_rshp ) )...
        .* bimag_sc_rshp;

    gradc = gradc_unsc / avgGlobalSARconstr;
end
end