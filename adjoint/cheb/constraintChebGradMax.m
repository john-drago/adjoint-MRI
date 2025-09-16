function [ c, gradc ] = constraintChebGradMax( pSc, opt )

gradMaxConstr = opt.gradMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

tdom = opt.tdom;

grad_idx = opt.grad_idx;
grad_sc = opt.scVec( grad_idx );
grad = grad_sc .* pSc( grad_idx );
grad_rshp = reshape( grad, [ opt.numCheb_grad, 3 ] );
[ grad_minmax, grad_minmax_pos ] = chebMinMax( grad_rshp, tdom );
[ grad_absmax, grad_absmax_rowpos ] = max( abs( grad_minmax ), [], 1 );

c_unsc = grad_absmax.'- gradMaxConstr;
c = c_unsc / gradMaxConstr;

if nargout > 1

    grad_absmax_linpos = sub2ind( size(grad_minmax), grad_absmax_rowpos, 1:3 );
    grad_absmax_pos = grad_minmax_pos( grad_absmax_linpos ).';

    T_n =  transpose(...
        evalChebClenshaw( grad_absmax_pos, eye( opt.numCheb_grad ), tdom ) );
    T_n_vec = reshape( T_n, [ opt.numCheb_grad * 3, 1 ] );
    
    gradc_unsc = zeros( opt.numVars, 3 );
    
    gradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        grad_idx, reshape( repmat( 1:3, [ opt.numCheb_grad, 1 ] ), [ opt.numCheb_grad * 3, 1 ] )  );

    grad_minmax_sign_vec = reshape(...
        repmat(...
        sign( grad_minmax( grad_absmax_linpos ) ), [ opt.numCheb_grad, 1 ] ),...
        [ opt.numCheb_grad * 3, 1 ] );

    gradc_unsc( gradc_linpos ) =...
        ( grad_minmax_sign_vec .* T_n_vec ) .* grad_sc;

    gradc = gradc_unsc / gradMaxConstr;

end

end