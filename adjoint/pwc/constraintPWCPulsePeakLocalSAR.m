function [ c, gradc ] = constraintPWCPulsePeakLocalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

% dutyCycle = opt.dutyCycle;
numXYCoils = opt.numXYCoils;
numTimePoints = opt.numTimePoints;
VOPs = opt.VOPs;

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

% maxTimeidx = 0;
% peakSAR = 0;
% highPowerVOP = 0;
% 
% for tt = 1:numTimePoints
% 
%     peakSARperVOP_tt = real( squeeze( sum( VOPs .* ( bcomp_rshp( tt, : )' * bcomp_rshp( tt, : ) ),...
%         [ 1 2 ] ) ) );
% 
%     [ peakSAR_tt, highPowerVOP_tt ] = max( peakSARperVOP_tt );
% 
%     if peakSAR_tt > peakSAR
%         peakSAR = peakSAR_tt;
%         maxTimeidx = tt;
%         highPowerVOP = highPowerVOP_tt;
%     end
% end

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

    % Determine variables indices for pTx channels corresponding to peak
    % timepoint
    var_timepoint_idx = transpose( maxTimeidx : numTimePoints : ( numXYCoils * numTimePoints ) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_timepoint_idx ) ) = real( gradf_zbar ) .* breal_sc( var_timepoint_idx );
    gradc_unsc( bimag_idx( var_timepoint_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_timepoint_idx );

    gradc = gradc_unsc / peakLocalSARconstr;
end
end