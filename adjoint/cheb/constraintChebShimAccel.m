function [ c, gradc ] = constraintChebShimAccel( pSc, opt )

shimAccelConstr = opt.shimAccel_constr;

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
ddshim_rshp = opt.D2_shim * shim_rshp;
[ ddshim_minmax, ddshim_minmax_pos ] = chebMinMax( ddshim_rshp, tdom );
[ ddshim_absmax, ddshim_absmax_rowpos ] = max( abs( ddshim_minmax ), [], 1 );

c_unsc = ddshim_absmax.'- shimAccelConstr;
c = c_unsc / shimAccelConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, numZCoils );

    ddshim_absmax_linpos = sub2ind( size(ddshim_minmax), ddshim_absmax_rowpos, 1:numZCoils );
    ddshim_absmax_pos = ddshim_minmax_pos( ddshim_absmax_linpos ).';

    T_n = transpose(...
         evalChebClenshaw( ddshim_absmax_pos, eye( opt.numCheb_shim - 2 ), tdom ) * opt.D2_shim );
    T_n_vec = reshape( T_n, [ ( opt.numCheb_shim ) * numZCoils, 1 ] );

    % ddshim_idx_rshp = shim_idx_rshp( 1:end, : );
    % ddshim_sc_rshp = shim_sc_rshp( 1:end, : );
    
    dgradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        shim_idx, reshape( repmat( 1:numZCoils, [ opt.numCheb_shim, 1 ] ), [ (opt.numCheb_shim) * numZCoils, 1 ] )  );

    dgrad_minmax_sign_vec = reshape(...
        repmat(...
        sign( ddshim_minmax( ddshim_absmax_linpos ) ), [ opt.numCheb_shim, 1 ] ),...
        [ ( opt.numCheb_shim ) * numZCoils, 1 ] );

    gradc_unsc( dgradc_linpos ) =...
        ( dgrad_minmax_sign_vec .* T_n_vec ) .* shim_sc;

    gradc = gradc_unsc / shimAccelConstr;

end

end