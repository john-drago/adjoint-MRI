function [ c, gradc ] = constraintPWCPulseAvgGlobalSAR( pSc, opt )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numTimePoints, numXYCoils ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numTimePoints, numXYCoils ] );

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
    gradf_zbar = reshape( gradf_zbar_rshp, [ numTimePoints * numXYCoils, 1 ] );


    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx ) = real( gradf_zbar ) .* breal_sc;
    gradc_unsc( bimag_idx ) = imag( gradf_zbar ) .* bimag_sc;


    gradc = gradc_unsc / avgGlobalSARconstr;
end
end