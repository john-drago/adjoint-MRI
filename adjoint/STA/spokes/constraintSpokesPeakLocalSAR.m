function [ c, gradc ] = constraintSpokesPeakLocalSAR( pSc, opt, ~ )
% This function will calculate the peak Local SAR for spoke pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

numSpokesCurr = length( opt.breal_idx ) / opt.numXYCoils;
numXYCoils = opt.numXYCoils;
VOPs = opt.VOPs;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, numSpokesCurr ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, numSpokesCurr ] );

bcomp_rshp = transpose( complex( breal_rshp, bimag_rshp ) );

if ~ismatrix( VOPs )

    localSARarr = real(...
        sum( conj( transpose( bcomp_rshp ) ) .* ...
        permute(...
        tensorprod( VOPs, transpose( bcomp_rshp ), 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) );

        if numSpokesCurr == 1
            localSARarr = squeeze( localSARarr );
        else
            localSARarr = permute( squeeze( localSARarr ), [ 2 1 ] );
        end
else

    localSARarr = transpose( real( ( VOPs * conj( bcomp_rshp ) ) .* bcomp_rshp ) );

end

[ peaklSAR, peaklSAR_li ] = max( localSARarr, [], 'all' );

c_unsc = ( peaklSAR ) - peakLocalSARconstr;
c = c_unsc / peakLocalSARconstr;

if nargout > 1

    % Determine peak VOP and time point
    [ highPowerVOP, maxSpokeidx ] = ind2sub( size(localSARarr), peaklSAR_li );
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    if ~ismatrix( VOPs )
        gradf_zbar = 2 * squeeze( VOPs( :, :, highPowerVOP ) ) * transpose( bcomp_rshp( maxSpokeidx, : ) );
    else
        gradf_zbar = 2 * VOPs * transpose( bcomp_rshp( maxSpokeidx, : ) );
    end

    % Determine variables indices for pTx channels corresponding to peak
    % spoke
    var_Spokes_idx = (maxSpokeidx-1) * numXYCoils + transpose( (1:numXYCoils) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_Spokes_idx ) ) = real( gradf_zbar ) .* breal_sc( var_Spokes_idx );
    gradc_unsc( bimag_idx( var_Spokes_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_Spokes_idx );

    gradc = gradc_unsc / peakLocalSARconstr;
end

end