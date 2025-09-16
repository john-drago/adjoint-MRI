function [ c, gradc ] = constraintPWCPulsePeakH10m( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakH10mconstr = opt.peakH10m_constr;

% dutyCycle = opt.dutyCycle;
numXYCoils = opt.numXYCoils;
numTimePoints = opt.numTimePoints;
VOPs_H10m = opt.VOPs_H10m;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numTimePoints, numXYCoils ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numTimePoints, numXYCoils ] );

bcomp_rshp_t = complex( breal_rshp, bimag_rshp );
bcomp_rshp = transpose( bcomp_rshp_t );

if ~ismatrix( VOPs_H10m )

    H10mArr = permute( sqrt( squeeze( real(...
        sum( conj( bcomp_rshp ) .* ...
        permute(...
        tensorprod( VOPs_H10m, bcomp_rshp, 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) ) ) ), [ 2 1 ] );

else

    H10mArr = transpose( sqrt( real( ( VOPs_H10m * conj( bcomp_rshp_t ) ) .* bcomp_rshp_t ) ) );

end

[ peakH10m, peakH10m_li ] = max( H10mArr, [], 'all' );

c_unsc = ( peakH10m ) - peakH10mconstr;
c = c_unsc / peakH10mconstr;

if nargout > 1

    % Determine peak VOP and time point
    [ highPowerVOP, maxTimeidx ] = ind2sub( size(H10mArr), peakH10m_li );
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    if all( peakH10m == 0 )
        peakH10m = 1;
    end
    if ~ismatrix( VOPs_H10m )
        gradf_zbar = ( 1/peakH10m ) * ( squeeze( VOPs_H10m( :, :, highPowerVOP ) ) * ( bcomp_rshp( :, maxTimeidx ) ) );
    else
        gradf_zbar = ( 1/peakH10m ) * ( VOPs_H10m * ( bcomp_rshp( :, maxTimeidx ) ) );
    end

    % Determine variables indices for pTx channels corresponding to peak
    % timepoint
    var_timepoint_idx = transpose( maxTimeidx : numTimePoints : ( numXYCoils * numTimePoints ) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_timepoint_idx ) ) = real( gradf_zbar ) .* breal_sc( var_timepoint_idx );
    gradc_unsc( bimag_idx( var_timepoint_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_timepoint_idx );

    gradc = gradc_unsc / peakH10mconstr;
end
end