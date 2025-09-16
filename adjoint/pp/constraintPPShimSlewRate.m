function [ c, gradc ] = constraintPPShimSlewRate( pSc, opt )

shimSlewRateConstr = opt.shimSlewRate_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numZCoils = opt.numZCoils;
PPStartStop_shim = opt.PPStartStop_shim;
numVarsPerChannel_shim = opt.numVarsPerChannel_shim;
varsToChebByPeriods_shim = opt.varsToChebByPeriods_shim;
numPPShape_shim = opt.numPPShape_shim;
numPPPeriods_shim = opt.numPPPeriods_shim;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
shim_rshp = reshape( shim, [ numVarsPerChannel_shim, numZCoils ] );
shimcheb_rshp = varsToChebByPeriods_shim * shim_rshp;
shimcheb_arr = permute( reshape( shimcheb_rshp, [ numPPShape_shim, numPPPeriods_shim, numZCoils ] ), [ 1, 3, 2 ] );
dshimcheb_arr = pagemtimes( opt.D_shim, shimcheb_arr );

dshim_minmax_arr = zeros( 2, numZCoils, numPPPeriods_shim );
dshim_minmax_pos_arr = zeros( 2, numZCoils, numPPPeriods_shim );
dshim_absmax_arr = zeros( numZCoils, numPPPeriods_shim );
dshim_absmax_rowpos_arr = zeros( numZCoils, numPPPeriods_shim );

for pp = 1:numPPPeriods_shim
    [ dshim_minmax_arr( :, :, pp ), dshim_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( dshimcheb_arr( :, :, pp ), PPStartStop_shim( pp, : ) );

    [ dshim_absmax_arr( :, pp ), dshim_absmax_rowpos_arr( :, pp ) ] = max( abs( dshim_minmax_arr( :, :, pp ) ), [], 1 );
end

dshim_minmax_arr = permute( dshim_minmax_arr, [ 1, 3, 2 ] );
dshim_minmax_pos_arr = permute( dshim_minmax_pos_arr, [ 1, 3, 2 ] );
dshim_absmax_rowpos_arr = transpose( dshim_absmax_rowpos_arr );
dshim_absmax_arr = transpose( dshim_absmax_arr );

dshim_absmax = reshape( dshim_absmax_arr, [ numPPPeriods_shim * numZCoils, 1 ] );

c_unsc = dshim_absmax - shimSlewRateConstr;
c = c_unsc / shimSlewRateConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, numZCoils * numPPPeriods_shim );

    numShim_numPPPeriods_arr = repmat( 1:numZCoils, [ numPPPeriods_shim, 1 ] );
    numPPPeriods_numShim_arr = repmat( transpose( 1:numPPPeriods_shim ), [ 1, numZCoils ] );

    dshim_absmax_linpos = sub2ind( size(dshim_minmax_arr ),...
        dshim_absmax_rowpos_arr, numPPPeriods_numShim_arr, numShim_numPPPeriods_arr );
    dshim_absmax_pos = dshim_minmax_pos_arr( dshim_absmax_linpos );
    
    for pp = 1:numPPPeriods_shim
        
        shim_minmax_Tn_samples = ...
            evalChebClenshaw( dshim_absmax_pos( pp, : ), eye( numPPShape_shim-1 ), PPStartStop_shim( pp, : ) );
        
        shimVarIdxsLength = opt.PPVarIdxs_shim(pp,2)-opt.PPVarIdxs_shim(pp,1)+1;
        shimVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_shim{ pp };
        shimVarCoilIdxShift = opt.shim_idx( shimVarCoilLocalIdxShift );
        constShimIdxShift = repmat( ( pp:numPPPeriods_shim:(numPPPeriods_shim*numZCoils) ) , [ shimVarIdxsLength, 1 ] );

        shapefunvals_shim = transpose( shim_minmax_Tn_samples *...
            ( opt.D_shim( :, :, pp ) * opt.shapeFnChebCoeffs_shim( :, opt.PPVarSFIdxs_shim(pp,1):opt.PPVarSFIdxs_shim(pp,2) ) ) );

        gradc_shim_linpos = sub2ind(...
            size( gradc_unsc ), shimVarCoilIdxShift, constShimIdxShift  );
        shim_minmax_sign_vec = repmat(...
            sign( dshim_minmax_arr( dshim_absmax_linpos(pp,:) ) ), [ shimVarIdxsLength, 1 ] );

        gradc_unsc( gradc_shim_linpos ) = ...
            ( shim_minmax_sign_vec .* shapefunvals_shim ) .* shim_sc( shimVarCoilLocalIdxShift );

    end    

    gradc = gradc_unsc / shimSlewRateConstr;

end

end