function [ c, gradc ] = constraintMPpTxPulseAvgGlobalSAR( pSc, opt )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
num_ORSP = opt.num_ORSP;
num_MPSP = opt.num_MPSP;
QGlobal = opt.QGlobal;

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );
brealORSP_rshp = reshape( brealORSP, [ num_ORSP, numXYCoils ] );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );
brealMPSP_rshp = reshape( brealMPSP, [ num_MPSP, numXYCoils ] );

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );
bimagORSP_rshp = reshape( bimagORSP, [ num_ORSP, numXYCoils ] );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );
bimagMPSP_rshp = reshape( bimagMPSP, [ num_MPSP, numXYCoils ] );

bcomp_rshp = complex(...
    [ brealORSP_rshp; brealMPSP_rshp ], [ bimagORSP_rshp; bimagMPSP_rshp ] );

b_pwr = bcomp_rshp' * bcomp_rshp;

avgGlobalSAR = ( (dutyCycle) / numTimePoints ) * real( sum( QGlobal .* b_pwr, [ 1 2 ] ) );

c_unsc = avgGlobalSAR - avgGlobalSARconstr;
% c = c_unsc / avgGlobalSARconstr;
c = c_unsc;

if nargout > 1

    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )

    % Determine contribution to average global SAR from each point
    gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
        * ( bcomp_rshp * conj( squeeze( QGlobal ) ) );
    
    gradf_zbar_ORSP = reshape( gradf_zbar_rshp( 1:num_ORSP, : ), [ num_ORSP * numXYCoils, 1 ] );
    gradf_zbar_MPSP = reshape( gradf_zbar_rshp( (num_ORSP+(1:num_MPSP)), : ), [ num_MPSP * numXYCoils, 1 ] );

    % Assign real and imaginary parts of the gradient
    gradc_unsc( brealORSP_idx ) = real( gradf_zbar_ORSP ) .* brealORSP_sc;
    gradc_unsc( bimagORSP_idx ) = imag( gradf_zbar_ORSP ) .* bimagORSP_sc;
    
    gradc_unsc( brealMPSP_idx ) = real( gradf_zbar_MPSP ) .* brealMPSP_sc;
    gradc_unsc( bimagMPSP_idx ) = imag( gradf_zbar_MPSP ) .* bimagMPSP_sc;

    % gradc = gradc_unsc / avgGlobalSARconstr;
    gradc = gradc_unsc;
end
end