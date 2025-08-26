function [ c, gradc ] = constraintChebShimMax( pSc, opt )

shimMaxConstr = opt.shimMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numZCoils = opt.numZCoils;
tdom = opt.tdom;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
shim_rshp = reshape( shim, [ opt.numCheb_shim, numZCoils ] );
[ shim_minmax, shim_minmax_pos ] = chebMinMax( shim_rshp, tdom );
[ shim_absmax, shim_absmax_rowpos ] = max( abs( shim_minmax ), [], 1 );

c_unsc = shim_absmax.'- shimMaxConstr;
c = c_unsc / shimMaxConstr;

if nargout > 1

    shim_absmax_linpos = sub2ind( size(shim_minmax), shim_absmax_rowpos, 1:numZCoils );
    shim_absmax_pos = shim_minmax_pos( shim_absmax_linpos ).';

    T_n =  transpose(...
        evalChebClenshaw( shim_absmax_pos, eye( opt.numCheb_shim ), tdom ) );
    T_n_vec = reshape( T_n, [ opt.numCheb_shim * numZCoils, 1 ] );
    
    gradc_unsc = zeros( opt.numVars, numZCoils );
    
    gradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        shim_idx, reshape( repmat( 1:numZCoils, [ opt.numCheb_shim, 1 ] ), [ opt.numCheb_shim * numZCoils, 1 ] )  );

    shim_minmax_sign_vec = reshape(...
        repmat(...
        sign( shim_minmax( shim_absmax_linpos ) ), [ opt.numCheb_shim, 1 ] ),...
        [ opt.numCheb_shim * numZCoils, 1 ] );

    gradc_unsc( gradc_linpos ) =...
        ( shim_minmax_sign_vec .* T_n_vec ) .* shim_sc;

    gradc = gradc_unsc / shimMaxConstr;

end

end