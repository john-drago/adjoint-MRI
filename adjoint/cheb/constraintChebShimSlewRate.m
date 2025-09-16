function [ c, gradc ] = constraintChebShimSlewRate( pSc, opt )

shimSlewRateConstr = opt.shimSlewRate_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

tdom = opt.tdom;
numZCoils = opt.numZCoils;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
% shim_idx_rshp = reshape( shim_idx, [ opt.numCheb_shim, numZCoils ]  );
% shim_sc_rshp = reshape( shim_sc, [ opt.numCheb_shim, numZCoils ] );
shim_rshp = reshape( shim, [ opt.numCheb_shim, numZCoils ] );
dshim_rshp = opt.D_shim * shim_rshp;
[ dshim_minmax, dshim_minmax_pos ] = chebMinMax( dshim_rshp, tdom );
[ dshim_absmax, dshim_absmax_rowpos ] = max( abs( dshim_minmax ), [], 1 );

c_unsc = dshim_absmax.'- shimSlewRateConstr;
c = c_unsc / shimSlewRateConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, numZCoils );

    dshim_absmax_linpos = sub2ind( size(dshim_minmax), dshim_absmax_rowpos, 1:numZCoils );
    dshim_absmax_pos = dshim_minmax_pos( dshim_absmax_linpos ).';

    T_n = transpose(...
         evalChebClenshaw( dshim_absmax_pos, eye( opt.numCheb_shim - 1 ), tdom ) * opt.D_shim );
    T_n_vec = reshape( T_n, [ ( opt.numCheb_shim ) * numZCoils, 1 ] );

    % dshim_idx_rshp = shim_idx_rshp( 1:end, : );
    % dshim_sc_rshp = shim_sc_rshp( 1:end, : );
    
    dgradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        shim_idx, reshape( repmat( 1:numZCoils, [ opt.numCheb_shim, 1 ] ), [ (opt.numCheb_shim) * numZCoils, 1 ] )  );

    dgrad_minmax_sign_vec = reshape(...
        repmat(...
        sign( dshim_minmax( dshim_absmax_linpos ) ), [ opt.numCheb_shim, 1 ] ),...
        [ ( opt.numCheb_shim ) * numZCoils, 1 ] );

    gradc_unsc( dgradc_linpos ) =...
        ( dgrad_minmax_sign_vec .* T_n_vec ) .* shim_sc;

    gradc = gradc_unsc / shimSlewRateConstr;

end

end