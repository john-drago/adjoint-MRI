function [ c, gradc ] = constraintChebPeakGlobalSAR( pSc, opt )
% This function will calculate the peak Global SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numCheb = opt.numCheb_RF;
TnT = opt.Tn( :, 1:numCheb );
QGlobal = opt.QGlobal;

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numCheb, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numCheb, numXYCoils ] );
breal_cheb = breal_sc .* pSc( breal_idx );
breal_cheb_rshp = reshape( breal_cheb, [ numCheb, numXYCoils ] );
breal_rshp = TnT * breal_cheb_rshp;

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numCheb, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numCheb, numXYCoils ] );
bimag_cheb = bimag_sc .* pSc( bimag_idx );
bimag_cheb_rshp = reshape( bimag_cheb, [ numCheb, numXYCoils ] );
bimag_rshp = TnT * bimag_cheb_rshp;

bcomp_rshp_t = complex( breal_rshp, bimag_rshp );
bcomp_rshp = transpose( bcomp_rshp_t );

globalSARarr = permute( squeeze( real(...
    sum( conj(bcomp_rshp) .* ... 
    tensorprod( QGlobal, bcomp_rshp, 2, 1 ) , 1 ) ) ), [ 2 1 ] );

[ peakGlobalSAR, maxTimeidx ] = max( globalSARarr, [], 'all' );

c_unsc = ( peakGlobalSAR ) - peakGlobalSARconstr;
c = c_unsc / peakGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    gradf_zbar = 2 * QGlobal * ( bcomp_rshp( :, maxTimeidx ) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx_rshp ) = transpose( real( gradf_zbar ) * TnT( maxTimeidx, :) ) .* breal_sc_rshp;
    gradc_unsc( bimag_idx_rshp ) = transpose( imag( gradf_zbar ) * TnT( maxTimeidx, :) ) .* bimag_sc_rshp;

    gradc = gradc_unsc / peakGlobalSARconstr;
end
end