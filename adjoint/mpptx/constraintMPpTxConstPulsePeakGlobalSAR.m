function [ c, gradc ] = constraintMPpTxConstPulsePeakGlobalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

% numXYCoils = opt.numXYCoils;
% numTimePoints = opt.numTimePoints;
% num_ORSP = opt.num_ORSP;
% num_MPSP = opt.num_MPSP;
QGlobal = opt.QGlobal;

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );

breal_rshp = [ brealORSP, brealMPSP ];

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );

bimag_rshp = [ bimagORSP, bimagMPSP ];

bcomp_rshp = complex( breal_rshp, bimag_rshp ).';

globalSARarr = permute( squeeze( real(...
    sum( conj( transpose( bcomp_rshp ) ) .* ... 
    tensorprod( QGlobal, transpose( bcomp_rshp ), 2, 1 ) , 1 ) ) ), [ 2 1 ] );

[ peakGlobalSAR, maxTimeidx ] = max( globalSARarr, [], 'all' );

c_unsc = ( peakGlobalSAR ) - peakGlobalSARconstr;
c = c_unsc / peakGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    gradf_zbar = 2 * QGlobal * ( bcomp_rshp( maxTimeidx, : ) ).';

    if maxTimeidx == 1
        breal_idx = brealORSP_idx;
        bimag_idx = bimagORSP_idx;
        breal_sc = brealORSP_sc;
        bimag_sc = bimagORSP_sc;
    elseif maxTimeidx == 2
        breal_idx = brealMPSP_idx;
        bimag_idx = bimagMPSP_idx;
        breal_sc = brealMPSP_sc;
        bimag_sc = bimagMPSP_sc;
    end
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx ) = real( gradf_zbar ) .* breal_sc;
    gradc_unsc( bimag_idx ) = imag( gradf_zbar ) .* bimag_sc;

    gradc = gradc_unsc / peakGlobalSARconstr;
end
end