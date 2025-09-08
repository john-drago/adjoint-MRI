function [ c, gradc ] = constraintkTPPulsePeakGlobalSAR( pSc, opt )
% This function will calculate the peak Global SAR for kT-point pulse and
% (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

peakGlobalSARconstr = opt.peakGlobalSAR_constr;

% dutyCycle = opt.dutyCycle;
% numTimePoints = opt.numTimePoints;
% RFSlewIntPtNum = opt.RFSlewIntPtNum;
% RFIntPtNum = opt.RFIntPtNum;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;

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
% peakGlobalSAR = 0;
% 
% for nn = 1:num_kTP
%     peakGlobalSAR_kTP = real( squeeze( sum( QGlobal .* ( bcomp_rshp( nn, : )' * bcomp_rshp( nn, : ) ),...
%         [ 1 2 ] ) ) );
% 
% 
%     if peakGlobalSAR_kTP > peakGlobalSAR
%         peakGlobalSAR = peakGlobalSAR_kTP;
%         maxkTidx = nn;
%     end
% end

globalSARarr = permute( squeeze( real(...
    sum( conj( transpose( bcomp_rshp ) ) .* ... 
    tensorprod( QGlobal, transpose( bcomp_rshp ), 2, 1 ) , 1 ) ) ), [ 2 1 ] );

[ peakGlobalSAR, maxkTidx ] = max( globalSARarr, [], 'all' );

c_unsc = ( peakGlobalSAR ) - peakGlobalSARconstr;
c = c_unsc / peakGlobalSARconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    gradf_zbar = 2 * QGlobal * transpose( ( bcomp_rshp( maxkTidx, : ) ) );

    % Determine variables indices for pTx channels corresponding to peak
    % kT-point
    var_kTP_idx = (maxkTidx-1) * numXYCoils + transpose( (1:numXYCoils) );
    
    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx( var_kTP_idx ) ) = real( gradf_zbar ) .* breal_sc( var_kTP_idx );
    gradc_unsc( bimag_idx( var_kTP_idx ) ) = imag( gradf_zbar ) .* bimag_sc( var_kTP_idx );

    gradc = gradc_unsc / peakGlobalSARconstr;
end
end