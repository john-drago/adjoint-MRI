function [ c, gradc ] = constraintMPpTxPulseAvgLocalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
num_ORSP = opt.num_ORSP;
num_MPSP = opt.num_MPSP;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_idx_rshp = reshape( brealORSP_idx, [ num_ORSP, numXYCoils ] );
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP_sc_rshp = reshape( brealORSP_sc, [ num_ORSP, numXYCoils ] );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );
brealORSP_rshp = reshape( brealORSP, [ num_ORSP, numXYCoils ] );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_idx_rshp = reshape( brealMPSP_idx, [ num_MPSP, numXYCoils ] );
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP_sc_rshp = reshape( brealMPSP_sc, [ num_MPSP, numXYCoils ] );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );
brealMPSP_rshp = reshape( brealMPSP, [ num_MPSP, numXYCoils ] );

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_idx_rshp = reshape( bimagORSP_idx, [ num_ORSP, numXYCoils ] );
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP_sc_rshp = reshape( bimagORSP_sc, [ num_ORSP, numXYCoils ] );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );
bimagORSP_rshp = reshape( bimagORSP, [ num_ORSP, numXYCoils ] );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_idx_rshp = reshape( bimagMPSP_idx, [ num_MPSP, numXYCoils ] );
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP_sc_rshp = reshape( bimagMPSP_sc, [ num_MPSP, numXYCoils ] );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );
bimagMPSP_rshp = reshape( bimagMPSP, [ num_MPSP, numXYCoils ] );

bcomp_rshp = complex(...
    [ brealORSP_rshp; brealMPSP_rshp ], [ bimagORSP_rshp; bimagMPSP_rshp ] );

b_pwr = ctranspose( bcomp_rshp ) * bcomp_rshp;

if ~ismatrix( VOPs )
    avgSAR = real( squeeze( sum( VOPs .* b_pwr, [ 1 2 ] ) ) );
else
    avgSAR = real( VOPs .* b_pwr );
end

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgSAR ) - avgLocalSARconstr;
c = c_unsc / avgLocalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, opt.numVOPs );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    if ~ismatrix( VOPs )
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) ) *...
            tensorprod( bcomp_rshp, conj( VOPs ), 2, 1 );
    else
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
            * ( bcomp_rshp * conj( VOPs ) );
    end

    gradf_zbar_rshp_ORSP = gradf_zbar_rshp( 1:num_ORSP, :, : );
    gradf_zbar_rshp_MPSP = gradf_zbar_rshp( (num_ORSP+1):end, :, : );

    gradf_zbar_brealORSP_rshp = real( gradf_zbar_rshp_ORSP ).* brealORSP_sc_rshp;
    gradf_zbar_bimagORSP_rshp = imag( gradf_zbar_rshp_ORSP ).* bimagORSP_sc_rshp;

    gradf_zbar_brealMPSP_rshp = real( gradf_zbar_rshp_MPSP ).* brealMPSP_sc_rshp;
    gradf_zbar_bimagMPSP_rshp = imag( gradf_zbar_rshp_MPSP ).* bimagMPSP_sc_rshp;

    % Assign real and imaginary parts of the gradient
    brealORSP_idx_rshp_rep = reshape( repmat( brealORSP_idx_rshp, [ 1, 1, numVOPs ] ), [ num_ORSP*numXYCoils, numVOPs ] );
    bimagORSP_idx_rshp_rep = reshape( repmat( bimagORSP_idx_rshp, [ 1, 1, numVOPs ] ), [ num_ORSP*numXYCoils, numVOPs ] );
    vop_bORSPidx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ num_ORSP, numXYCoils ] ), [ num_ORSP*numXYCoils, numVOPs ] );
    lin_brealORSP_idx = sub2ind( size( gradc_unsc ), brealORSP_idx_rshp_rep(:), vop_bORSPidx_rshp_rep(:) );
    lin_bimagORSP_idx = sub2ind( size( gradc_unsc ), bimagORSP_idx_rshp_rep(:), vop_bORSPidx_rshp_rep(:) );
    gradc_unsc( lin_brealORSP_idx ) = gradf_zbar_brealORSP_rshp( : );
    gradc_unsc( lin_bimagORSP_idx ) = gradf_zbar_bimagORSP_rshp( : );

    brealMPSP_idx_rshp_rep = reshape( repmat( brealMPSP_idx_rshp, [ 1, 1, numVOPs ] ), [ num_MPSP*numXYCoils, numVOPs ] );
    bimagMPSP_idx_rshp_rep = reshape( repmat( bimagMPSP_idx_rshp, [ 1, 1, numVOPs ] ), [ num_MPSP*numXYCoils, numVOPs ] );
    vop_bMPSPidx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ num_MPSP, numXYCoils ] ), [ num_MPSP*numXYCoils, numVOPs ] );
    lin_brealMPSP_idx = sub2ind( size( gradc_unsc ), brealMPSP_idx_rshp_rep(:), vop_bMPSPidx_rshp_rep(:) );
    lin_bimagMPSP_idx = sub2ind( size( gradc_unsc ), bimagMPSP_idx_rshp_rep(:), vop_bMPSPidx_rshp_rep(:) );
    gradc_unsc( lin_brealMPSP_idx ) = gradf_zbar_brealMPSP_rshp( : );
    gradc_unsc( lin_bimagMPSP_idx ) = gradf_zbar_bimagMPSP_rshp( : );

    gradc = gradc_unsc / avgLocalSARconstr;
end
end