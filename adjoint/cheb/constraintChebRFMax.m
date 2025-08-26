function [ c, gradc ] = constraintChebRFMax( pSc, opt )

RFMaxConstr = opt.RFMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numXYCoils = opt.numXYCoils;
tdom = opt.tdom;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ opt.numCheb_RF, numXYCoils ] );
[ breal_minmax, breal_minmax_pos ] = chebMinMax( breal_rshp, tdom );
[ breal_absmax, breal_absmax_rowpos ] = max( abs( breal_minmax ), [], 1 );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ opt.numCheb_RF, numXYCoils ] );
[ bimag_minmax, bimag_minmax_pos ] = chebMinMax( bimag_rshp, tdom );
[ bimag_absmax, bimag_absmax_rowpos ] = max( abs( bimag_minmax ), [], 1 );

c_unsc = [ transpose( breal_absmax ); transpose( bimag_absmax ) ] - RFMaxConstr;
c = c_unsc / RFMaxConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2*numXYCoils );

    breal_absmax_linpos = sub2ind( size(breal_minmax), breal_absmax_rowpos, 1:numXYCoils );
    breal_absmax_pos = transpose( breal_minmax_pos( breal_absmax_linpos ) );
    
    T_n_breal =  transpose(...
        evalChebClenshaw( breal_absmax_pos, eye( opt.numCheb_RF ), tdom ) );
    T_n_breal_vec = reshape( T_n_breal, [ opt.numCheb_RF * numXYCoils, 1 ] );

    gradc_breal_linpos = sub2ind(...
        size( gradc_unsc ),...
        breal_idx, reshape( repmat( 1:numXYCoils, [ opt.numCheb_RF, 1 ] ), [ opt.numCheb_RF * numXYCoils, 1 ] )  );
    breal_minmax_sign_vec = reshape(...
        repmat(...
        sign( breal_minmax( breal_absmax_linpos ) ), [ opt.numCheb_RF, 1 ] ),...
        [ opt.numCheb_RF * numXYCoils, 1 ] );

    bimag_absmax_linpos = sub2ind( size(bimag_minmax), bimag_absmax_rowpos, 1:numXYCoils );
    bimag_absmax_pos = transpose( bimag_minmax_pos( bimag_absmax_linpos ) );
    
    T_n_bimag =  transpose(...
        evalChebClenshaw( bimag_absmax_pos, eye( opt.numCheb_RF ), tdom ) );
    T_n_bimag_vec = reshape( T_n_bimag, [ opt.numCheb_RF * numXYCoils, 1 ] );

    gradc_bimag_linpos = sub2ind(...
        size( gradc_unsc ),...
        bimag_idx, reshape( repmat( numXYCoils + (1:numXYCoils), [ opt.numCheb_RF, 1 ] ), [ opt.numCheb_RF * numXYCoils, 1 ] )  );
    bimag_minmax_sign_vec = reshape(...
        repmat(...
        sign( bimag_minmax( bimag_absmax_linpos ) ), [ opt.numCheb_RF, 1 ] ),...
        [ opt.numCheb_RF * numXYCoils, 1 ] );

    gradc_unsc( gradc_breal_linpos ) =...
        ( breal_minmax_sign_vec .* T_n_breal_vec ) .* breal_sc;
    gradc_unsc( gradc_bimag_linpos ) =...
        ( bimag_minmax_sign_vec .* T_n_bimag_vec ) .* bimag_sc;

    gradc = gradc_unsc / RFMaxConstr;

end

end