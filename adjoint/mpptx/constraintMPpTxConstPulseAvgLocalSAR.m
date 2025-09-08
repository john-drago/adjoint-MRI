function [ c, gradc ] = constraintMPpTxConstPulseAvgLocalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
% num_ORSP = opt.num_ORSP;
% num_MPSP = opt.num_MPSP;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;
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

bcomp_rshp = transpose( complex( breal_rshp, bimag_rshp ) );

b_pwr_unscale = bcomp_rshp' * bcomp_rshp;
b_const_pwr = bcomp_rshp' * ( bcomp_rshp .* [ tORSPIntPtNum; tMPSPIntPtNum ] );
b_slew_pwr = sum( b_pwr_unscale .* reshape( slewScale.^2, [ 1, 1, RFSlewIntPtNum ] ) , 3 );

b_pwr = 2 * b_slew_pwr + b_const_pwr;

if ~ismatrix( VOPs )
    avgLocalSAR = real( squeeze( sum( VOPs .* b_pwr, [ 1 2 ] ) ) );
else
    avgLocalSAR = real( VOPs .* b_pwr );
end

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgLocalSAR ) - avgLocalSARconstr;
c = c_unsc / avgLocalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numVOPs );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    if ~ismatrix( VOPs )

        % gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 * ( [ tORSPIntPtNum, tMPSPIntPtNum ] + 2 * sum( slewScale.^2 ) ) ) )...
        %     .* ( squeeze( VOPs( :, :, peakVOPlocation ) ) * ( bcomp_rshp ).' );

        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) ) *...
            tensorprod( ( bcomp_rshp .* ( [ tORSPIntPtNum; tMPSPIntPtNum ] + 2 * sum( slewScale.^2 ) ) ), conj( VOPs ), 2, 1 );
        gradf_zbar_ORSP = squeeze( gradf_zbar_rshp( 1, :, : ) );
        gradf_zbar_MPSP = squeeze( gradf_zbar_rshp( 2, :, : ) );

    else
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 * ( [ tORSPIntPtNum, tMPSPIntPtNum ] + 2 * sum( slewScale.^2 ) ) ) )...
            .* (  VOPs * ( bcomp_rshp ).' );
        gradf_zbar_ORSP = gradf_zbar_rshp( 1 );
        gradf_zbar_MPSP = gradf_zbar_rshp( 2 );
    end

    gradf_zbar_brealORSP = real( gradf_zbar_ORSP ).* brealORSP_sc;
    gradf_zbar_bimagORSP = imag( gradf_zbar_ORSP ).* bimagORSP_sc;

    gradf_zbar_brealMPSP = real( gradf_zbar_MPSP ).* brealMPSP_sc;
    gradf_zbar_bimagMPSP = imag( gradf_zbar_MPSP ).* bimagMPSP_sc;

    % Assign real and imaginary parts of the gradient
    brealORSP_idx_rep = repmat( brealORSP_idx, [ 1, numVOPs ] );
    bimagORSP_idx_rep = repmat( bimagORSP_idx, [ 1, numVOPs ] );
    vop_bORSPidx_rep = repmat( reshape( 1:numVOPs, [ 1, numVOPs  ] ), [ numXYCoils, 1 ] );
    lin_brealORSP_idx = sub2ind( size( gradc_unsc ), brealORSP_idx_rep(:), vop_bORSPidx_rep(:) );
    lin_bimagORSP_idx = sub2ind( size( gradc_unsc ), bimagORSP_idx_rep(:), vop_bORSPidx_rep(:) );
    gradc_unsc( lin_brealORSP_idx ) = gradf_zbar_brealORSP( : );
    gradc_unsc( lin_bimagORSP_idx ) = gradf_zbar_bimagORSP( : );
    
    brealMPSP_idx_rep = repmat( brealMPSP_idx, [ 1, numVOPs ] );
    bimagMPSP_idx_rep = repmat( bimagMPSP_idx, [ 1, numVOPs ] );
    vop_bMPSPidx_rep = repmat( reshape( 1:numVOPs, [ 1, numVOPs  ] ), [ numXYCoils, 1 ] );
    lin_brealMPSP_idx = sub2ind( size( gradc_unsc ), brealMPSP_idx_rep(:), vop_bMPSPidx_rep(:) );
    lin_bimagMPSP_idx = sub2ind( size( gradc_unsc ), bimagMPSP_idx_rep(:), vop_bMPSPidx_rep(:) );
    gradc_unsc( lin_brealMPSP_idx ) = gradf_zbar_brealMPSP( : );
    gradc_unsc( lin_bimagMPSP_idx ) = gradf_zbar_bimagMPSP( : );

    gradc = gradc_unsc / avgLocalSARconstr;
end
end