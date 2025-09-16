function [ c, gradc ] = constraintChebAvgLocalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numCheb = opt.numCheb_RF;
TnT = opt.Tn( :, 1:numCheb );
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;

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

        gradf_zbar_breal_rshp = reshape( tensorprod( TnT, real( gradf_zbar_rshp ), 1, 1) .* breal_sc_rshp, [ numCheb*numXYCoils, numVOPs ] );
        gradf_zbar_bimag_rshp = reshape( tensorprod( TnT, imag( gradf_zbar_rshp ), 1, 1) .* bimag_sc_rshp, [ numCheb*numXYCoils, numVOPs ] );

    else
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * ( 2 ) )...
            * ( bcomp_rshp * conj( VOPs ) );

        gradf_zbar_breal_rshp = ( transpose( TnT ) * real( gradf_zbar_rshp ) ) .* breal_sc_rshp;
        gradf_zbar_bimag_rshp = ( transpose( TnT ) * imag( gradf_zbar_rshp ) ) .* bimag_sc_rshp;
    end

    % Assign real and imaginary parts of the gradient
    breal_idx_rshp_rep = reshape( repmat( breal_idx_rshp, [ 1, 1, numVOPs ] ), [ numCheb*numXYCoils, numVOPs ] );
    bimag_idx_rshp_rep = reshape( repmat( bimag_idx_rshp, [ 1, 1, numVOPs ] ), [ numCheb*numXYCoils, numVOPs ] );
    vop_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ numCheb, numXYCoils ] ), [ numCheb*numXYCoils, numVOPs ] );
    lin_breal_idx = sub2ind( size( gradc_unsc ), breal_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    lin_bimag_idx = sub2ind( size( gradc_unsc ), bimag_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    gradc_unsc( lin_breal_idx ) = gradf_zbar_breal_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradf_zbar_bimag_rshp( : );

    gradc = gradc_unsc / avgLocalSARconstr;
end
end