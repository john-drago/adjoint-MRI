function [ c, gradc ] = constraintFourierGradSlewRate( pSc, opt, maxExtrema )

if nargin < 3
    maxExtrema = 1;
else
    maxExtrema = min( [ maxExtrema, opt.orderFourier_grad - 1 ] );
    maxExtrema = max( [ 1, maxExtrema ] );
end

gradSlewRateConstr = opt.gradSlewRate_constr;

% Assume pSc is scaled to be between -1 and 1
pSc = pSc( : );

numCoils = 3;
numFourier = opt.numFourier_grad;
FBgradMat = opt.D_grad;

grad_fourier_idx = opt.grad_idx;
grad_fourier_idx_rshp = reshape( grad_fourier_idx, [ numFourier, numCoils ] );
grad_fourier_sc = opt.scVec( grad_fourier_idx );
grad_fourier_sc_rshp = reshape( grad_fourier_sc, [ numFourier, numCoils ] );
grad_fourier = grad_fourier_sc .* pSc( grad_fourier_idx );
grad_fourier_rshp = reshape( grad_fourier, [ numFourier, numCoils ] );

grad = FBgradMat * grad_fourier_rshp;
grad_abs = abs( grad );

[ grad_max, grad_max_pos ] = maxk( grad_abs, maxExtrema, 1 );

c_unsc = reshape( grad_max, [ maxExtrema*3, 1 ] ) - gradSlewRateConstr;
c = c_unsc / gradSlewRateConstr;

if nargout > 1

    numCoilsRepMat = repmat( 1:numCoils, [ maxExtrema, 1 ] );

    grad_max_vals_idx = sub2ind( size( grad_abs ),...
        grad_max_pos, numCoilsRepMat );
    grad_max_vals_idx_rshp = reshape( grad_max_vals_idx, [ maxExtrema*numCoils, 1 ] );
    grad_max_pos_rshp = reshape( grad_max_pos, [ maxExtrema*numCoils, 1 ] );
    % numCoilsRepMat_rshp = reshape( numCoilsRepMat, [ maxExtrema*numCoils, 1 ] );

    gradc_unsc = zeros( opt.numVars, numCoils * maxExtrema );

    dgradabsdgrad = sign( grad( grad_max_vals_idx_rshp ) );
    dgradabsdgradf = dgradabsdgrad .* FBgradMat( grad_max_pos_rshp, : );

    grad_fourier_idx_subpos = reshape( permute( repmat( grad_fourier_idx_rshp, [ 1, 1, maxExtrema ] ), [ 1, 3, 2 ] ), [ numFourier, maxExtrema*numCoils ] );
    grad_coil_idx_subpos = repmat( uint32( 1:(numCoils*maxExtrema) ), [ numFourier, 1 ] );

    gradc_gradabs_gradf_linpos = sub2ind( size( gradc_unsc ),...
        grad_fourier_idx_subpos,...
        grad_coil_idx_subpos );

    gradc_unsc( gradc_gradabs_gradf_linpos ) =...
        transpose( dgradabsdgradf ) .* reshape( permute( repmat( grad_fourier_sc_rshp, [ 1, 1, maxExtrema ] ), [ 1, 3, 2 ] ), [ numFourier, maxExtrema*numCoils ] );

    gradc = gradc_unsc / gradSlewRateConstr;


end

end