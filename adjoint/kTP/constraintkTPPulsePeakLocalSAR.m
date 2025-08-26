function [ c, gradc ] = constraintkTPPulsePeakLocalSAR( pSc, opt )
% This function will calculate the peak Local SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakLocalSARconstr = opt.peakLocalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
% RFSlewIntPtNum = opt.RFSlewIntPtNum;
% RFIntPtNum = opt.RFIntPtNum;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
VOPs = opt.VOPs;

% slewScale = ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, num_kTP ] );
% breal_slew = repmat( breal, [ 1, RFSlewIntPtNum ] ) .* ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );
% breal_slew_rshp = permute( reshape( breal_slew, [ numXYCoils, num_kTP, RFSlewIntPtNum ] ), [ 2 1 3 ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
% bimag_slew = repmat( bimag, [ 1, RFSlewIntPtNum ] ) .* ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );
% bimag_slew_rshp = permute( reshape( bimag_slew, [ numXYCoils, num_kTP, RFSlewIntPtNum ] ), [ 2 1 3 ] );

bcomp_rshp = transpose( complex( breal_rshp, bimag_rshp ) );
% bcomp_slew_rshp = breal_slew_rshp + 1j * bimag_slew_rshp;

% maxkTidx = 0;
% peakSAR = 0;
% highPowerVOP = 0;
% 
% for nn = 1:num_kTP
% 
%     peakSARperVOP_kTP = real( squeeze( sum( VOPs .* ( bcomp_rshp( nn, : )' * bcomp_rshp( nn, : ) ),...
%         [ 1 2 ] ) ) );
% 
%     % bcomp_rshp_ktp = repmat( bcomp_rshp( :, nn ), [ 1 1 size( VOPs, 3 ) ] );
%     % peakSARperVOP_kTP = real( pagemtimes( pagemtimes( pagectranspose( bcomp_rshp_ktp ), VOPs ),  bcomp_rshp_ktp ) );
% 
%     [ peakSAR_kTP, highPowerVOP_kTP ] = max( peakSARperVOP_kTP );
% 
%     if peakSAR_kTP > peakSAR
%         peakSAR = peakSAR_kTP;
%         maxkTidx = nn;
%         highPowerVOP = highPowerVOP_kTP;
%     end
% end

if ~ismatrix( VOPs )

    localSARarr = squeeze( real(...
        sum( conj( transpose( bcomp_rshp ) ) .* ...
        permute(...
        tensorprod( VOPs, transpose( bcomp_rshp ), 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) ) );

    if num_kTP > 1
        localSARarr = transpose( localSARarr );
    end
    
else

    localSARarr = transpose( real( ( VOPs * conj( bcomp_rshp ) ) .* bcomp_rshp ) );

end

[ peaklSAR, peaklSAR_li ] = max( localSARarr, [], 'all' );

c_unsc = ( peaklSAR ) - peakLocalSARconstr;
c = c_unsc / peakLocalSARconstr;

if nargout > 1

    % Determine peak VOP and time point
    [ highPowerVOP, maxkTidx ] = ind2sub( size(localSARarr), peaklSAR_li );
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    if ~ismatrix( VOPs )
        gradf_zbar = 2 * squeeze( VOPs( :, :, highPowerVOP ) ) * transpose( bcomp_rshp( maxkTidx, : ) );
    else
        gradf_zbar = 2 * VOPs * transpose( bcomp_rshp( maxkTidx, : ) );
    end

    % Determine variables indices for pTx channels corresponding to peak
    % kT-point
    var_kTP_idx = (maxkTidx-1) * numXYCoils + (1:numXYCoils).';
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_kTP_idx ) ) = real( gradf_zbar ) .* breal_sc( var_kTP_idx );
    gradc_unsc( bimag_idx( var_kTP_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_kTP_idx );

    gradc = gradc_unsc / peakLocalSARconstr;
end
end