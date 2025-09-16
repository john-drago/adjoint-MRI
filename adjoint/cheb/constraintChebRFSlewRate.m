function [ c, gradc ] = constraintChebRFSlewRate( pSc, opt )

RFSlewRateConstr = opt.RFSlewRate_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numXYCoils = opt.numXYCoils;
tdom = opt.tdom;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
% breal_idx_rshp = reshape( breal_idx, [ opt.numCheb_RF, numXYCoils ]  );
% breal_sc_rshp = reshape( breal_sc, [ opt.numCheb_RF, numXYCoils ] );
breal_rshp = reshape( breal, [ opt.numCheb_RF, numXYCoils ] );
dbreal_rshp = opt.D_RF * breal_rshp;
[ dbreal_minmax, dbreal_minmax_pos ] = chebMinMax( dbreal_rshp, tdom );
[ dbreal_absmax, dbreal_absmax_rowpos ] = max( abs( dbreal_minmax ), [], 1 );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
% bimag_idx_rshp = reshape( bimag_idx, [ opt.numCheb_RF, numXYCoils ]  );
% bimag_sc_rshp = reshape( bimag_sc, [ opt.numCheb_RF, numXYCoils ] );
bimag_rshp = reshape( bimag, [ opt.numCheb_RF, numXYCoils ] );
dbimag_rshp = opt.D_RF * bimag_rshp;
[ dbimag_minmax, dbimag_minmax_pos ] = chebMinMax( dbimag_rshp, tdom );
[ dbimag_absmax, dbimag_absmax_rowpos ] = max( abs( dbimag_minmax ), [], 1 );

c_unsc = [ dbreal_absmax.'; dbimag_absmax.' ] - RFSlewRateConstr;
c = c_unsc / RFSlewRateConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2*numXYCoils );

    dbreal_absmax_linpos = sub2ind( size(dbreal_minmax), dbreal_absmax_rowpos, 1:numXYCoils );
    dbreal_absmax_pos = dbreal_minmax_pos( dbreal_absmax_linpos ).';
    
    T_n_dbreal =  transpose(...
        evalChebClenshaw( dbreal_absmax_pos, eye( opt.numCheb_RF-1 ), tdom ) * opt.D_RF );
    T_n_dbreal_vec = reshape( T_n_dbreal, [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    % dbreal_idx_rshp = breal_idx_rshp( 1:end, : );
    % dbreal_sc_rshp = breal_sc_rshp( 1:end, : );

    gradc_dbreal_linpos = sub2ind(...
        size( gradc_unsc ),...
        breal_idx, reshape( repmat( 1:numXYCoils, [ ( opt.numCheb_RF ) , 1 ] ), [ ( opt.numCheb_RF )  * numXYCoils, 1 ] )  );
    dbreal_minmax_sign_vec = reshape(...
        repmat(...
        sign( dbreal_minmax( dbreal_absmax_linpos ) ), [ ( opt.numCheb_RF ), 1 ] ),...
        [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    dbimag_absmax_linpos = sub2ind( size(dbimag_minmax), dbimag_absmax_rowpos, 1:numXYCoils );
    dbimag_absmax_pos = dbimag_minmax_pos( dbimag_absmax_linpos ).';
    
    T_n_dbimag =  transpose(...
        evalChebClenshaw( dbimag_absmax_pos, eye( opt.numCheb_RF-1 ), tdom ) * opt.D_RF );
    T_n_dbimag_vec = reshape( T_n_dbimag, [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    % dbimag_idx_rshp = bimag_idx_rshp( 1:end, : );
    % dbimag_sc_rshp = bimag_sc_rshp( 1:end, : );

    gradc_dbimag_linpos = sub2ind(...
        size( gradc_unsc ),...
        bimag_idx, reshape( repmat( numXYCoils + (1:numXYCoils), [ ( opt.numCheb_RF ) , 1 ] ), [ ( opt.numCheb_RF ) * numXYCoils, 1 ] )  );
    dbimag_minmax_sign_vec = reshape(...
        repmat(...
        sign( dbimag_minmax( dbimag_absmax_linpos ) ), [ ( opt.numCheb_RF ) , 1 ] ),...
        [ ( opt.numCheb_RF ) * numXYCoils, 1 ] );

    gradc_unsc( gradc_dbreal_linpos ) =...
        ( dbreal_minmax_sign_vec .* T_n_dbreal_vec ) .* breal_sc;
    gradc_unsc( gradc_dbimag_linpos ) =...
        ( dbimag_minmax_sign_vec .* T_n_dbimag_vec ) .* bimag_sc;

    gradc = gradc_unsc / RFSlewRateConstr;

end

end