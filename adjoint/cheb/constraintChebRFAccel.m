function [ c, gradc ] = constraintChebRFAccel( pSc, opt )

RFAccelConstr = opt.RFAccel_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

tdom = opt.tdom;
numXYCoils = opt.numXYCoils;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
% breal_idx_rshp = reshape( breal_idx, [ opt.numCheb_RF, numXYCoils ]  );
% breal_sc_rshp = reshape( breal_sc, [ opt.numCheb_RF, numXYCoils ] );
breal_rshp = reshape( breal, [ opt.numCheb_RF, numXYCoils ] );
ddbreal_rshp = opt.D2_RF * breal_rshp;
[ ddbreal_minmax, ddbreal_minmax_pos ] = chebMinMax( ddbreal_rshp, tdom );
[ ddbreal_absmax, ddbreal_absmax_rowpos ] = max( abs( ddbreal_minmax ), [], 1 );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
% bimag_idx_rshp = reshape( bimag_idx, [ opt.numCheb_RF, numXYCoils ]  );
% bimag_sc_rshp = reshape( bimag_sc, [ opt.numCheb_RF, numXYCoils ] );
bimag_rshp = reshape( bimag, [ opt.numCheb_RF, numXYCoils ] );
ddbimag_rshp = opt.D2_RF * bimag_rshp;
[ ddbimag_minmax, ddbimag_minmax_pos ] = chebMinMax( ddbimag_rshp, tdom );
[ ddbimag_absmax, ddbimag_absmax_rowpos ] = max( abs( ddbimag_minmax ), [], 1 );

c_unsc = [ ddbreal_absmax.'; ddbimag_absmax.' ] - RFAccelConstr;
c = c_unsc / RFAccelConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2*numXYCoils );

    ddbreal_absmax_linpos = sub2ind( size(ddbreal_minmax), ddbreal_absmax_rowpos, 1:numXYCoils );
    ddbreal_absmax_pos = ddbreal_minmax_pos( ddbreal_absmax_linpos ).';

    T_n_breal = transpose(...
         evalChebClenshaw( ddbreal_absmax_pos, eye( opt.numCheb_RF - 2 ), tdom ) * opt.D2_RF );
    T_n_breal_vec = reshape( T_n_breal, [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    % ddbreal_idx_rshp = breal_idx_rshp( 1:end, : );
    % ddbreal_sc_rshp = breal_sc_rshp( 1:end, : );
    
    dgradc_breal_linpos = sub2ind(...
        size( gradc_unsc ),...
        breal_idx, reshape( repmat( 1:numXYCoils, [ opt.numCheb_RF, 1 ] ), [ (opt.numCheb_RF) * numXYCoils, 1 ] )  );

    dgrad_breal_minmax_sign_vec = reshape(...
        repmat(...
        sign( ddbreal_minmax( ddbreal_absmax_linpos ) ), [ opt.numCheb_RF, 1 ] ),...
        [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    ddbimag_absmax_linpos = sub2ind( size(ddbimag_minmax), ddbimag_absmax_rowpos, 1:numXYCoils );
    ddbimag_absmax_pos = ddbimag_minmax_pos( ddbimag_absmax_linpos ).';

    T_n_bimag = transpose(...
         evalChebClenshaw( ddbimag_absmax_pos, eye( opt.numCheb_RF - 2 ), tdom ) * opt.D2_RF );
    T_n_bimag_vec = reshape( T_n_bimag, [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    % ddbimag_idx_rshp = bimag_idx_rshp( 1:end, : );
    % ddbimag_sc_rshp = bimag_sc_rshp( 1:end, : );
    
    dgradc_bimag_linpos = sub2ind(...
        size( gradc_unsc ),...
        bimag_idx, reshape( repmat( numXYCoils + (1:numXYCoils), [ opt.numCheb_RF, 1 ] ), [ (opt.numCheb_RF) * numXYCoils, 1 ] )  );

    dgrad_bimag_minmax_sign_vec = reshape(...
        repmat(...
        sign( ddbimag_minmax( ddbimag_absmax_linpos ) ), [ opt.numCheb_RF, 1 ] ),...
        [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    gradc_unsc( dgradc_breal_linpos ) =...
        ( dgrad_breal_minmax_sign_vec .* T_n_breal_vec ) .* breal_sc;
    gradc_unsc( dgradc_bimag_linpos ) =...
        ( dgrad_bimag_minmax_sign_vec .* T_n_bimag_vec ) .* bimag_sc;

    gradc = gradc_unsc / RFAccelConstr;

end

end