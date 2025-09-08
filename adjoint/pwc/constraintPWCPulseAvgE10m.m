function [ c, gradc ] = constraintPWCPulseAvgE10m( pSc, opt )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgE10mconstr = opt.avgE10m_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
VOPs_E10m = opt.VOPs_E10m;
numVOPs_E10m = uint32( opt.numVOPs_E10m );

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numTimePoints, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numTimePoints, numXYCoils ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numTimePoints, numXYCoils ] );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numTimePoints, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numTimePoints, numXYCoils ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numTimePoints, numXYCoils ] );

bcomp_rshp = complex( breal_rshp, bimag_rshp );

b_pwr = ctranspose( bcomp_rshp ) * bcomp_rshp;

if ~ismatrix( VOPs_E10m )
    avgE10m = sqrt( real( squeeze( sum( VOPs_E10m .* b_pwr, [ 1 2 ] ) ) ) );
else
    avgE10m = sqrt( real( VOPs_E10m .* b_pwr ) );
end

c_unsc = ( (dutyCycle) / numTimePoints ) * ( avgE10m ) - avgE10mconstr;
c = c_unsc / avgE10mconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numVOPs_E10m );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    if all( avgE10m == 0, "all" )
        avgE10m = ones( numVOPs_E10m, 1 );
    end
    if ~ismatrix( VOPs_E10m )
        gradf_zbar_rshp = ( ( (dutyCycle) / numTimePoints ) * reshape( 1./avgE10m, [ 1, 1, numVOPs_E10m ] ) ) .*...
            tensorprod( bcomp_rshp, conj( VOPs_E10m ), 2, 1 );
    else
        gradf_zbar_rshp = ( ( (dutyCycle) / ( numTimePoints * avgE10m ) ) )...
            * ( bcomp_rshp * conj( VOPs_E10m ) );
    end

    gradf_zbar_breal_rshp = real( gradf_zbar_rshp ).* breal_sc_rshp;
    gradf_zbar_bimag_rshp = imag( gradf_zbar_rshp ).* bimag_sc_rshp;

    % Assign real and imaginary parts of the gradient
    breal_idx_rshp_rep = reshape( repmat( breal_idx_rshp, [ 1, 1, numVOPs_E10m ] ),...
        [ numTimePoints*numXYCoils, numVOPs_E10m ] );
    bimag_idx_rshp_rep = reshape( repmat( bimag_idx_rshp, [ 1, 1, numVOPs_E10m ] ),...
        [ numTimePoints*numXYCoils, numVOPs_E10m ] );
    vop_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs_E10m, [ 1, 1, numVOPs_E10m  ] ),...
        [ numTimePoints, numXYCoils ] ), [ numTimePoints*numXYCoils, numVOPs_E10m ] );
    lin_breal_idx = sub2ind( size( gradc_unsc ),...
        breal_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    lin_bimag_idx = sub2ind( size( gradc_unsc ),...
        bimag_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
    gradc_unsc( lin_breal_idx ) = gradf_zbar_breal_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradf_zbar_bimag_rshp( : );

    gradc = gradc_unsc / avgE10mconstr;
end
end