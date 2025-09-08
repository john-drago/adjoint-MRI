function [ c, gradc ] = constraintMPpTxPulsePeakLocalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

numXYCoils = opt.numXYCoils;
% numTimePoints = opt.numTimePoints;
num_ORSP = opt.num_ORSP;
num_MPSP = opt.num_MPSP;
VOPs = opt.VOPs;

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );
brealORSP_rshp = reshape( brealORSP, [ num_ORSP, numXYCoils ] );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );
brealMPSP_rshp = reshape( brealMPSP, [ num_MPSP, numXYCoils ] );

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );
bimagORSP_rshp = reshape( bimagORSP, [ num_ORSP, numXYCoils ] );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );
bimagMPSP_rshp = reshape( bimagMPSP, [ num_MPSP, numXYCoils ] );

bcomp_rshp_t = complex(...
    [ brealORSP_rshp; brealMPSP_rshp ], [ bimagORSP_rshp; bimagMPSP_rshp ] );
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
    if maxTimeidx > num_ORSP
        ORSPMPSP_maxidx = maxTimeidx - num_ORSP;
        % Determine variables indices for pTx channels corresponding to peak
        % timepoint
        var_timepoint_idx = ( ORSPMPSP_maxidx : num_MPSP : ( numXYCoils * num_MPSP ) ).';
        breal_idx = brealMPSP_idx( var_timepoint_idx );
        bimag_idx = bimagMPSP_idx( var_timepoint_idx );
        breal_sc = brealMPSP_sc( var_timepoint_idx );
        bimag_sc = bimagMPSP_sc( var_timepoint_idx );
    else
        ORSPMPSP_maxidx = maxTimeidx;
        % Determine variables indices for pTx channels corresponding to peak
        % timepoint
        var_timepoint_idx = ( ORSPMPSP_maxidx : num_ORSP : ( numXYCoils * num_ORSP ) ).';
        breal_idx = brealORSP_idx( var_timepoint_idx );
        bimag_idx = bimagORSP_idx( var_timepoint_idx );
        breal_sc = brealORSP_sc( var_timepoint_idx );
        bimag_sc = bimagORSP_sc( var_timepoint_idx );
    end
    
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
    gradc_unsc( breal_idx ) = real( gradf_zbar ) .* breal_sc;
    gradc_unsc( bimag_idx ) = imag( gradf_zbar ) .* bimag_sc;

    gradc = gradc_unsc / peakLocalSARconstr;
end
end