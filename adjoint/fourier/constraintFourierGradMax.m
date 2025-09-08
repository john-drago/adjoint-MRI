function [ c, gradc ] = constraintFourierGradMax( pSc, opt )

gradMaxConstr = opt.gradMax_constr;

% Assume pSc is scaled to be between -1 and 1
pSc = pSc( : );

numCoils = 3;
numFourier = opt.numFourier_grad;
FBgradMat = opt.FBgrad;

grad_fourier_idx = opt.grad_idx;
grad_fourier_idx_rshp = reshape( grad_fourier_idx, [ numFourier, numCoils ] );
grad_fourier_sc = opt.scVec( grad_fourier_idx );
grad_fourier_sc_rshp = reshape( grad_fourier_sc, [ numFourier, numCoils ] );
grad_fourier = grad_fourier_sc .* pSc( grad_fourier_idx );
grad_fourier_rshp = reshape( grad_fourier, [ numFourier, numCoils ] );

grad = FBgradMat * grad_fourier_rshp;
grad_abs = abs( grad );

[ grad_max, grad_max_pos ] = max( grad_abs, [], 1 );

c_unsc = transpose( grad_max ) - gradMaxConstr;
c = c_unsc / gradMaxConstr;

if nargout > 1

    grad_max_vals_idx = sub2ind( size( grad_abs ),...
        grad_max_pos, 1:numCoils );

    gradc_unsc = zeros( opt.numVars, numCoils );

    dgradabsdgrad = sign( grad( grad_max_vals_idx ) );
    dgradabsdgradf = transpose( dgradabsdgrad ) .* FBgradMat( grad_max_pos, : );
    
    gradc_gradabs_gradf = sub2ind( size( gradc_unsc ),...
        transpose( grad_fourier_idx_rshp ),...
        repmat( transpose(1:numCoils), [ 1, numFourier ] ) );
    gradc_unsc( gradc_gradabs_gradf ) =...
        ( dgradabsdgradf ) .* transpose( grad_fourier_sc_rshp );

    gradc = gradc_unsc / gradMaxConstr;

end

end