function [ c, gradc ] = constraintChebGradAccel( pSc, opt )

gradAccelConstr = opt.gradAccel_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

tdom = opt.tdom;
% numZCoils = opt.numZCoils;

grad_idx = opt.grad_idx;
grad_sc = opt.scVec( grad_idx );
grad = grad_sc .* pSc( grad_idx );
% grad_idx_rshp = reshape( grad_idx, [ opt.numCheb_grad, numZCoils ]  );
% grad_sc_rshp = reshape( grad_sc, [ opt.numCheb_grad, numZCoils ] );
grad_rshp = reshape( grad, [ opt.numCheb_grad, 3 ] );
ddgrad_rshp = opt.D2_grad * grad_rshp;
[ ddgrad_minmax, ddgrad_minmax_pos ] = chebMinMax( ddgrad_rshp, tdom );
[ ddgrad_absmax, ddgrad_absmax_rowpos ] = max( abs( ddgrad_minmax ), [], 1 );

c_unsc = ddgrad_absmax.'- gradAccelConstr;
c = c_unsc / gradAccelConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 3 );

    ddgrad_absmax_linpos = sub2ind( size(ddgrad_minmax), ddgrad_absmax_rowpos, 1:3 );
    ddgrad_absmax_pos = ddgrad_minmax_pos( ddgrad_absmax_linpos ).';

    T_n = transpose(...
         evalChebClenshaw( ddgrad_absmax_pos, eye( opt.numCheb_grad - 2 ), tdom ) * opt.D2_grad );
    T_n_vec = reshape( T_n, [ ( opt.numCheb_grad ) * 3, 1 ] );

    % ddgrad_idx_rshp = grad_idx_rshp( 1:end, : );
    % ddgrad_sc_rshp = grad_sc_rshp( 1:end, : );
    
    dgradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        grad_idx, reshape( repmat( 1:3, [ opt.numCheb_grad, 1 ] ), [ (opt.numCheb_grad) * 3, 1 ] )  );

    dgrad_minmax_sign_vec = reshape(...
        repmat(...
        sign( ddgrad_minmax( ddgrad_absmax_linpos ) ), [ opt.numCheb_grad, 1 ] ),...
        [ ( opt.numCheb_grad ) * 3, 1 ] );

    gradc_unsc( dgradc_linpos ) =...
        ( dgrad_minmax_sign_vec .* T_n_vec ) .* grad_sc;

    gradc = gradc_unsc / gradAccelConstr;

end

end