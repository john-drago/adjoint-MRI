function [ c, gradc ] = constraintChebRFBandwidth( pSc, opt )

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

RFBandwidthconstr = opt.RFBandwidth_constr;

adjDFTmtxTn = opt.adjDFTmtxTn;
% numAdjDFTmtx = size( adjDFTmtxTn, 1 );

numXYCoils = opt.numXYCoils;
numCheb = opt.numCheb_RF;

breal_idx = opt.breal_idx;
breal_idx_rshp = reshape( breal_idx, [ numCheb, numXYCoils ] );
breal_sc = opt.scVec( breal_idx );
breal_sc_rshp = reshape( breal_sc, [ numCheb, numXYCoils ] );
breal_cheb = breal_sc .* pSc( breal_idx );
breal_cheb_rshp = reshape( breal_cheb, [ numCheb, numXYCoils ] );

bimag_idx = opt.bimag_idx;
bimag_idx_rshp = reshape( bimag_idx, [ numCheb, numXYCoils ] );
bimag_sc = opt.scVec( bimag_idx );
bimag_sc_rshp = reshape( bimag_sc, [ numCheb, numXYCoils ] );
bimag_cheb = bimag_sc .* pSc( bimag_idx );
bimag_cheb_rshp = reshape( bimag_cheb, [ numCheb, numXYCoils ] );

bcomp_rshp = complex( breal_cheb_rshp, bimag_cheb_rshp );

bcomp_rshp_fft = adjDFTmtxTn * bcomp_rshp;

bcomp_rshp_fft_abs = abs( bcomp_rshp_fft );

c_unsc = transpose( opt.adjdf * sum( bcomp_rshp_fft_abs, 1 )  - RFBandwidthconstr );

c = c_unsc / RFBandwidthconstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numXYCoils );
    
    gradf_zbar_rshp = ...
         opt.adjdf * ctranspose( adjDFTmtxTn ) * ( bcomp_rshp_fft ./ bcomp_rshp_fft_abs ); 
    
    lin_breal_idx = sub2ind( [ opt.numVars, numXYCoils ],...
         breal_idx_rshp, repmat( 1:numXYCoils, [ numCheb, 1 ] ) );
    lin_bimag_idx = sub2ind( [ opt.numVars, numXYCoils ],...
         bimag_idx_rshp, repmat( 1:numXYCoils, [ numCheb, 1 ] ) );

    gradc_unsc( lin_breal_idx ) = real( gradf_zbar_rshp ) .* breal_sc_rshp;
    gradc_unsc( lin_bimag_idx ) = imag( gradf_zbar_rshp ) .* bimag_sc_rshp;

    gradc = gradc_unsc / RFBandwidthconstr;


end
end