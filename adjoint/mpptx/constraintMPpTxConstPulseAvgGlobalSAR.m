function [ c, gradc ] = constraintMPpTxConstPulseAvgGlobalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
% numXYCoils = opt.numXYCoils;
% num_ORSP = opt.num_ORSP;
% num_MPSP = opt.num_MPSP;
QGlobal = opt.QGlobal;
tORSPIntPtNum = opt.tORSPIntPtNum;
tMPSPIntPtNum = opt.tMPSPIntPtNum;
RFSlewIntPtNum = opt.RFSlewIntPtNum;

slewScale = ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );

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

b_pwr_unscale = bcomp_rshp' * bcomp_rshp;
b_const_pwr = bcomp_rshp' * ( bcomp_rshp .* [ tORSPIntPtNum; tMPSPIntPtNum ] );
b_slew_pwr = sum( b_pwr_unscale .* reshape( slewScale.^2, [ 1, 1, RFSlewIntPtNum ] ) , 3 );

b_pwr = 2 * b_slew_pwr + b_const_pwr;

avgGlobalSAR = real( squeeze( sum( QGlobal .* b_pwr, [ 1 2 ] ) ) );

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgGlobalSAR ) - avgGlobalSARconstr;
c = c_unsc / avgGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 * ( [ tORSPIntPtNum, tMPSPIntPtNum ] + 2 * sum( slewScale.^2 ) ) ) )...
        .* (  QGlobal * ( bcomp_rshp ).' );

    gradf_zbar_ORSP = gradf_zbar_rshp( :, 1 );
    gradf_zbar_MPSP = gradf_zbar_rshp( :, 2 );

    % Assign real and imaginary parts of the gradient
    gradc_unsc( brealORSP_idx ) = real( gradf_zbar_ORSP ) .* brealORSP_sc;
    gradc_unsc( bimagORSP_idx ) = imag( gradf_zbar_ORSP ) .* bimagORSP_sc;
    
    gradc_unsc( brealMPSP_idx ) = real( gradf_zbar_MPSP ) .* brealMPSP_sc;
    gradc_unsc( bimagMPSP_idx ) = imag( gradf_zbar_MPSP ) .* bimagMPSP_sc;

    gradc = gradc_unsc / avgGlobalSARconstr;
end
end