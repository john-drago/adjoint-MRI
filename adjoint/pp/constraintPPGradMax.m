function [ c, gradc ] = constraintPPGradMax( pSc, opt )

gradMaxConstr = opt.gradMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

PPStartStop_grad = opt.PPStartStop_grad;
numVarsPerChannel_grad = opt.numVarsPerChannel_grad;
varsToChebByPeriods_grad = opt.varsToChebByPeriods_grad;
numPPShape_grad = opt.numPPShape_grad;
numPPPeriods_grad = opt.numPPPeriods_grad;

grad_idx = opt.grad_idx;
grad_sc = opt.scVec( grad_idx );
grad = grad_sc .* pSc( grad_idx );
grad_rshp = reshape( grad, [ numVarsPerChannel_grad, 3 ] );
gradcheb_rshp = varsToChebByPeriods_grad * grad_rshp;
gradcheb_arr = permute( reshape( gradcheb_rshp, [ numPPShape_grad, numPPPeriods_grad, 3 ] ), [ 1, 3, 2 ] );

grad_minmax_arr = zeros( 2, 3, numPPPeriods_grad );
grad_minmax_pos_arr = zeros( 2, 3, numPPPeriods_grad );
grad_absmax_arr = zeros( 3, numPPPeriods_grad );
grad_absmax_rowpos_arr = zeros( 3, numPPPeriods_grad );

for pp = 1:numPPPeriods_grad
    [ grad_minmax_arr( :, :, pp ), grad_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( gradcheb_arr( :, :, pp ), PPStartStop_grad( pp, : ) );

    [ grad_absmax_arr( :, pp ), grad_absmax_rowpos_arr( :, pp ) ] = max( abs( grad_minmax_arr( :, :, pp ) ), [], 1 );
end

grad_minmax_arr = permute( grad_minmax_arr, [ 1, 3, 2 ] );
grad_minmax_pos_arr = permute( grad_minmax_pos_arr, [ 1, 3, 2 ] );
grad_absmax_rowpos_arr = transpose( grad_absmax_rowpos_arr );
grad_absmax_arr = transpose( grad_absmax_arr );

grad_absmax = reshape( grad_absmax_arr, [ numPPPeriods_grad * 3, 1 ] );

c_unsc = grad_absmax - gradMaxConstr;
c = c_unsc / gradMaxConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 3 * numPPPeriods_grad );

    numGrad_numPPPeriods_arr = repmat( 1:3, [ numPPPeriods_grad, 1 ] );
    numPPPeriods_numGrad_arr = repmat( transpose( 1:numPPPeriods_grad ), [ 1, 3 ] );

    grad_absmax_linpos = sub2ind( size(grad_minmax_arr ),...
        grad_absmax_rowpos_arr, numPPPeriods_numGrad_arr, numGrad_numPPPeriods_arr );
    grad_absmax_pos = grad_minmax_pos_arr( grad_absmax_linpos );
    
    for pp = 1:numPPPeriods_grad
        
        grad_minmax_Tn_samples = ...
            evalChebClenshaw( grad_absmax_pos( pp, : ), eye( numPPShape_grad ), PPStartStop_grad( pp, : ) );
        
        gradVarIdxsLength = opt.PPVarIdxs_grad(pp,2)-opt.PPVarIdxs_grad(pp,1)+1;
        gradVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_grad{ pp };
        gradVarCoilIdxShift = opt.grad_idx( gradVarCoilLocalIdxShift );
        constGradIdxShift = repmat( ( pp:numPPPeriods_grad:(numPPPeriods_grad*3) ) , [ gradVarIdxsLength, 1 ] );

        shapefunvals_grad = transpose( grad_minmax_Tn_samples * opt.shapeFnChebCoeffs_grad( :, opt.PPVarSFIdxs_grad(pp,1):opt.PPVarSFIdxs_grad(pp,2) )  );

        gradc_grad_linpos = sub2ind(...
            size( gradc_unsc ), gradVarCoilIdxShift, constGradIdxShift  );
        grad_minmax_sign_vec = repmat(...
            sign( grad_minmax_arr( grad_absmax_linpos(pp,:) ) ), [ gradVarIdxsLength, 1 ] );

        gradc_unsc( gradc_grad_linpos ) = ...
            ( grad_minmax_sign_vec .* shapefunvals_grad ) .* grad_sc( gradVarCoilLocalIdxShift );

    end    

    gradc = gradc_unsc / gradMaxConstr;

end

end