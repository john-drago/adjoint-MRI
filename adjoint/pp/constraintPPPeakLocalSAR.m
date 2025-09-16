function [ c, gradc ] = constraintPPPeakLocalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
varsToTimepoints_RF = opt.varsToTimepoints_RF;
VOPs = opt.VOPs;

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

bcomp_rshp_t = complex( breal_rshp, bimag_rshp );
bcomp_rshp = transpose( bcomp_rshp_t );

if ~ismatrix( VOPs )

    localSARarr = permute( squeeze( real(...
        sum( conj(bcomp_rshp) .* ...
        permute(...
        tensorprod( VOPs, bcomp_rshp, 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) ) ), [ 2 1 ] );

else

    localSARarr = transpose( real( ( VOPs * conj( bcomp_rshp_t ) ) .* bcomp_rshp_t ) );

end

[ peaklSAR, peaklSAR_li ] = max( localSARarr, [], 'all' );

c_unsc = ( peaklSAR ) - peakLocalSARconstr;
c = c_unsc / peakLocalSARconstr;

if nargout > 1

    % Determine peak VOP and time point
    [ highPowerVOP, maxTimeidx ] = ind2sub( size(localSARarr), peaklSAR_li );
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    if ~ismatrix( VOPs )
        gradf_zbar = 2 * squeeze( VOPs( :, :, highPowerVOP ) ) * ( bcomp_rshp( :, maxTimeidx ) );
    else
        gradf_zbar = 2 * VOPs * ( bcomp_rshp( :, maxTimeidx ) );
    end
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx_rshp ) = transpose( real( gradf_zbar ) * varsToTimepoints_RF( maxTimeidx, :) ) .* breal_sc_rshp;
    gradc_unsc( bimag_idx_rshp ) = transpose( imag( gradf_zbar ) * varsToTimepoints_RF( maxTimeidx, :) ) .* bimag_sc_rshp;

    gradc = gradc_unsc / peakLocalSARconstr;
end
end