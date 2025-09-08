function [ c, gradc ] = constraintFourierAvgLocalSAR( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numFourier = opt.numFourier_RF;
FBRFMat = opt.FBRF;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;

breal_fourier_idx = opt.breal_idx;
breal_fourier_idx_rshp = reshape( breal_fourier_idx, [ numFourier, numXYCoils ] );
breal_fourier_sc = opt.scVec( breal_fourier_idx );
breal_fourier_sc_rshp = reshape( breal_fourier_sc, [ numFourier, numXYCoils ] );
breal_fourier = breal_fourier_sc .* pSc( breal_fourier_idx );
breal_fourier_rshp = reshape( breal_fourier, [ numFourier, numXYCoils ] );

bimag_fourier_idx = opt.bimag_idx;
bimag_fourier_idx_rshp = reshape( bimag_fourier_idx, [ numFourier, numXYCoils ] );
bimag_fourier_sc = opt.scVec( bimag_fourier_idx );
bimag_fourier_sc_rshp = reshape( bimag_fourier_sc, [ numFourier, numXYCoils ] );
bimag_fourier = bimag_fourier_sc .* pSc( bimag_fourier_idx );
bimag_fourier_rshp = reshape( bimag_fourier, [ numFourier, numXYCoils ] );

bcomp_rshp = complex( breal_fourier_rshp, bimag_fourier_rshp );
bcomp = ( FBRFMat * bcomp_rshp );

b_pwr = ctranspose( bcomp ) * bcomp;

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
    grad_breal_brealf = real( FBRFMat );
    grad_breal_bimagf = -imag( FBRFMat );
    grad_bimag_brealf = imag( FBRFMat );
    grad_bimag_bimagf = real( FBRFMat );
    
    % Determine contribution to average local SAR from each point
    if ~ismatrix( VOPs )

        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) ) *...
            tensorprod( bcomp, conj( VOPs ), 2, 1 );

        gradf_zbar_breal_brealf_rshp =...
            reshape( tensorprod( grad_breal_brealf, real( gradf_zbar_rshp ), 1, 1) .* breal_fourier_sc_rshp, [ numFourier*numXYCoils, numVOPs ] );
        gradf_zbar_breal_bimagf_rshp =...
            reshape( tensorprod( grad_breal_bimagf, real( gradf_zbar_rshp ), 1, 1) .* bimag_fourier_sc_rshp, [ numFourier*numXYCoils, numVOPs ] );

        gradf_zbar_bimag_brealf_rshp =...
            reshape( tensorprod( grad_bimag_brealf, imag( gradf_zbar_rshp ), 1, 1) .* breal_fourier_sc_rshp, [ numFourier*numXYCoils, numVOPs ] );
        gradf_zbar_bimag_bimagf_rshp =...
            reshape( tensorprod( grad_bimag_bimagf, imag( gradf_zbar_rshp ), 1, 1) .* bimag_fourier_sc_rshp, [ numFourier*numXYCoils, numVOPs ] );
    else
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
            * ( bcomp * conj( VOPs ) );

        gradf_zbar_breal_brealf_rshp =... 
            ( transpose( grad_breal_brealf ) * real( gradf_zbar_rshp ) ) .* breal_fourier_sc_rshp;
        gradf_zbar_breal_bimagf_rshp =...
            ( transpose( grad_breal_bimagf ) * real( gradf_zbar_rshp ) ) .* bimag_fourier_sc_rshp;

        gradf_zbar_bimag_brealf_rshp =...
            ( transpose( grad_bimag_brealf ) * imag( gradf_zbar_rshp ) ) .* breal_fourier_sc_rshp;
        gradf_zbar_bimag_bimagf_rshp =...
            ( transpose( grad_bimag_bimagf ) * imag( gradf_zbar_rshp ) ) .* bimag_fourier_sc_rshp;
    end

    % Assign real and imaginary parts of the gradient
    breal_idx_rshp_rep = reshape( repmat( breal_fourier_idx_rshp, [ 1, 1, numVOPs ] ), [ numFourier*numXYCoils, numVOPs ] );
    bimag_idx_rshp_rep = reshape( repmat( bimag_fourier_idx_rshp, [ 1, 1, numVOPs ] ), [ numFourier*numXYCoils, numVOPs ] );
    vop_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ numFourier, numXYCoils ] ), [ numFourier*numXYCoils, numVOPs ] );
    lin_breal_idx = sub2ind( size( gradc_unsc ), breal_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    lin_bimag_idx = sub2ind( size( gradc_unsc ), bimag_idx_rshp_rep(:), vop_idx_rshp_rep(:) );

    gradc_unsc( lin_breal_idx ) = gradc_unsc( lin_breal_idx ) +...
        gradf_zbar_breal_brealf_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradc_unsc( lin_bimag_idx ) +...
        gradf_zbar_breal_bimagf_rshp( : );

    gradc_unsc( lin_breal_idx ) = gradc_unsc( lin_breal_idx ) +...
        gradf_zbar_bimag_brealf_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradc_unsc( lin_bimag_idx ) +...
        gradf_zbar_bimag_bimagf_rshp( : );

    gradc = gradc_unsc / avgLocalSARconstr;
    
end
end