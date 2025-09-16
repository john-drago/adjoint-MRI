function [ c, gradc ] = constraintSpokesPeakGlobalSAR( pSc, opt, ~ )
% This function will calculate the peak Global SAR for spokes pulses and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

numSpokesCurr = length( opt.breal_idx ) / opt.numXYCoils;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, numSpokesCurr ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, numSpokesCurr ] );

bcomp_rshp = transpose( complex( breal_rshp, bimag_rshp ) );

globalSARarr = permute( squeeze( real(...
    sum( conj( transpose( bcomp_rshp ) ) .* ... 
    tensorprod( QGlobal, transpose( bcomp_rshp ), 2, 1 ) , 1 ) ) ), [ 2 1 ] );

[ peakGlobalSAR, maxSpokeidx ] = max( globalSARarr, [], 'all' );

c_unsc = ( peakGlobalSAR ) - peakGlobalSARconstr;
c = c_unsc / peakGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    gradf_zbar = 2 * QGlobal * ( transpose( bcomp_rshp( maxSpokeidx, : ) ) );

    % Determine variables indices for pTx channels corresponding to peak
    % spoke
    var_spoke_idx = (maxSpokeidx-1) * numXYCoils + transpose( (1:numXYCoils) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_spoke_idx ) ) = real( gradf_zbar ) .* breal_sc( var_spoke_idx );
    gradc_unsc( bimag_idx( var_spoke_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_spoke_idx );

    gradc = gradc_unsc / peakGlobalSARconstr;
end

end