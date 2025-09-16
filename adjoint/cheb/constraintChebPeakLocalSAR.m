function [ c, gradc ] = constraintChebPeakLocalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
numXYCoils = opt.numXYCoils;
numCheb = opt.numCheb_RF;
TnT = opt.Tn( :, 1:numCheb );
VOPs = opt.VOPs;

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
    gradc_unsc( breal_idx_rshp ) = transpose( real( gradf_zbar ) * TnT( maxTimeidx, :) ) .* breal_sc_rshp;
    gradc_unsc( bimag_idx_rshp ) = transpose( imag( gradf_zbar ) * TnT( maxTimeidx, :) ) .* bimag_sc_rshp;

    gradc = gradc_unsc / peakLocalSARconstr;
end
end