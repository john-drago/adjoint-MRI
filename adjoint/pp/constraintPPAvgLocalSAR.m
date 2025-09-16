function [ c, gradc ] = constraintPPAvgLocalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
varsToTimepoints_RF = opt.varsToTimepoints_RF;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;

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
bimag_shape = bimag_sc .* pSc( bimag_idx );
bimag_shape_rshp = reshape( bimag_shape, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_rshp = varsToTimepoints_RF * bimag_shape_rshp;

bcomp_rshp = complex( breal_rshp, bimag_rshp );

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
        
        % gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) ) *...
        %     tensorprod( TnT, tensorprod( bcomp_rshp, conj( VOPs ), 2, 1 ), 1, 1 );

        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) ) *...
            tensorprod( bcomp_rshp, conj( VOPs ), 2, 1 );

        gradf_zbar_breal_rshp = reshape( tensorprod( varsToTimepoints_RF, real( gradf_zbar_rshp ), 1, 1) .* breal_sc_rshp, [ numVarsPerChannel_RF*numXYCoils, numVOPs ] );
        gradf_zbar_bimag_rshp = reshape( tensorprod( varsToTimepoints_RF, imag( gradf_zbar_rshp ), 1, 1) .* bimag_sc_rshp, [ numVarsPerChannel_RF*numXYCoils, numVOPs ] );
    else
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
            * ( bcomp_rshp * conj( VOPs ) );

        gradf_zbar_breal_rshp = ( transpose( varsToTimepoints_RF ) * real( gradf_zbar_rshp ) ) .* breal_sc_rshp;
        gradf_zbar_bimag_rshp = ( transpose( varsToTimepoints_RF ) * imag( gradf_zbar_rshp ) ) .* bimag_sc_rshp;
    end

    % Assign real and imaginary parts of the gradient
    breal_idx_rshp_rep = reshape( repmat( breal_idx_rshp, [ 1, 1, numVOPs ] ), [ numVarsPerChannel_RF*numXYCoils, numVOPs ] );
    bimag_idx_rshp_rep = reshape( repmat( bimag_idx_rshp, [ 1, 1, numVOPs ] ), [ numVarsPerChannel_RF*numXYCoils, numVOPs ] );
    vop_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ numVarsPerChannel_RF, numXYCoils ] ), [ numVarsPerChannel_RF*numXYCoils, numVOPs ] );
    lin_breal_idx = sub2ind( size( gradc_unsc ), breal_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    lin_bimag_idx = sub2ind( size( gradc_unsc ), bimag_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    gradc_unsc( lin_breal_idx ) = gradf_zbar_breal_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradf_zbar_bimag_rshp( : );

    gradc = gradc_unsc / avgLocalSARconstr;
end
end