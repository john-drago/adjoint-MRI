function [ c, gradc ] = constraintPWCPulsePeakGlobalSAR( pSc, opt )
% This function will calculate the peak Global SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

% dutyCycle = opt.dutyCycle;
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

bcomp_rshp_t = complex( breal_rshp, bimag_rshp );
bcomp_rshp = transpose( bcomp_rshp_t );

% maxTimeidx = 0;
% peakGlobalSAR = 0;
% 
% for nn = 1:num_kTP
%     peakGlobalSAR_tt = real( squeeze( sum( QGlobal .* ( bcomp_rshp( nn, : )' * bcomp_rshp( nn, : ) ),...
%         [ 1 2 ] ) ) );
% 
%     if peakGlobalSAR_tt > peakGlobalSAR
%         peakGlobalSAR = peakGlobalSAR_tt;
%         maxTimeidx = nn;
%     end
% end

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

    % Determine variables indices for pTx channels corresponding to peak
    % timepoint
    var_timepoint_idx = transpose( maxTimeidx : numTimePoints : ( numXYCoils * numTimePoints ) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_timepoint_idx ) ) = real( gradf_zbar ) .* breal_sc( var_timepoint_idx );
    gradc_unsc( bimag_idx( var_timepoint_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_timepoint_idx );

    gradc = gradc_unsc / peakGlobalSARconstr;
end
end