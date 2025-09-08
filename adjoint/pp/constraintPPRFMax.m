function [ c, gradc ] = constraintPPRFMax( pSc, opt )

RFMaxConstr = opt.RFMax_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numXYCoils = opt.numXYCoils;
PPStartStop_RF = opt.PPStartStop_RF;
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
varsToChebByPeriods_RF = opt.varsToChebByPeriods_RF;
numPPShape_RF = opt.numPPShape_RF;
numPPPeriods_RF = opt.numPPPeriods_RF;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numVarsPerChannel_RF, numXYCoils ] );
brealcheb_rshp = varsToChebByPeriods_RF * breal_rshp;
brealcheb_arr = permute( reshape( brealcheb_rshp, [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] ), [ 1, 3, 2 ] );

breal_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
breal_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
breal_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
breal_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numVarsPerChannel_RF, numXYCoils ] );
bimagcheb_rshp = varsToChebByPeriods_RF * bimag_rshp;
bimagcheb_arr = permute( reshape( bimagcheb_rshp, [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] ), [ 1, 3, 2 ] );

bimag_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
bimag_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
bimag_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
bimag_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

for pp = 1:numPPPeriods_RF
    [ breal_minmax_arr( :, :, pp ), breal_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( brealcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );
    [ bimag_minmax_arr( :, :, pp ), bimag_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( bimagcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );

    [ breal_absmax_arr( :, pp ), breal_absmax_rowpos_arr( :, pp ) ] = max( abs( breal_minmax_arr( :, :, pp ) ), [], 1 );
    [ bimag_absmax_arr( :, pp ), bimag_absmax_rowpos_arr( :, pp ) ] = max( abs( bimag_minmax_arr( :, :, pp ) ), [], 1 );
end

breal_minmax_arr = permute( breal_minmax_arr, [ 1, 3, 2 ] );
bimag_minmax_arr = permute( bimag_minmax_arr, [ 1, 3, 2 ] );

breal_minmax_pos_arr = permute( breal_minmax_pos_arr, [ 1, 3, 2 ] );
bimag_minmax_pos_arr = permute( bimag_minmax_pos_arr, [ 1, 3, 2 ] );

breal_absmax_rowpos_arr = transpose( breal_absmax_rowpos_arr );
bimag_absmax_rowpos_arr = transpose( bimag_absmax_rowpos_arr );

breal_absmax_arr = transpose( breal_absmax_arr );
bimag_absmax_arr = transpose( bimag_absmax_arr );

breal_absmax = reshape( breal_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );
bimag_absmax = reshape( bimag_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );

c_unsc = [ breal_absmax; bimag_absmax ] - RFMaxConstr;
c = c_unsc / RFMaxConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2 * numXYCoils * numPPPeriods_RF );

    numXYCoils_numPPPeriods_arr = repmat(  1:numXYCoils, [ numPPPeriods_RF, 1 ] );
    numPPPeriods_numXYCoils_arr = repmat( transpose( 1:numPPPeriods_RF ), [ 1, numXYCoils ] );

    % breal
    breal_absmax_linpos = sub2ind( size(breal_minmax_arr ),...
        breal_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    breal_absmax_pos = breal_minmax_pos_arr( breal_absmax_linpos );

    % bimag
    bimag_absmax_linpos = sub2ind( size(bimag_minmax_arr ),...
        bimag_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    bimag_absmax_pos = bimag_minmax_pos_arr( bimag_absmax_linpos );
    
    for pp = 1:numPPPeriods_RF
        
        % breal
        breal_minmax_Tn_samples = ...
            evalChebClenshaw( breal_absmax_pos( pp, : ), eye( numPPShape_RF ), PPStartStop_RF( pp, : ) );
        
        bVarIdxsLength = opt.PPVarIdxs_RF(pp,2)-opt.PPVarIdxs_RF(pp,1)+1;
        brealVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        brealVarCoilIdxShift = opt.breal_idx( brealVarCoilLocalIdxShift );
        constCoilIdxShift = repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );

        shapefunvals_breal = transpose( breal_minmax_Tn_samples * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) )  );

        gradc_breal_linpos = sub2ind(...
            size( gradc_unsc ), brealVarCoilIdxShift, constCoilIdxShift  );
        breal_minmax_sign_vec = repmat(...
            sign( breal_minmax_arr( breal_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_breal_linpos ) = ...
            ( breal_minmax_sign_vec .* shapefunvals_breal ) .* breal_sc( brealVarCoilLocalIdxShift );

        % bimag
        bimag_minmax_Tn_samples = ...
            evalChebClenshaw( bimag_absmax_pos( pp, : ), eye( numPPShape_RF ), PPStartStop_RF( pp, : ) );
        
        bimagVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        bimagVarCoilIdxShift = opt.bimag_idx( bimagVarCoilLocalIdxShift );
        constCoilIdxShift = numPPPeriods_RF * numXYCoils +...
            repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );
        
        shapefunvals_bimag = transpose( bimag_minmax_Tn_samples * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) )  );

        gradc_bimag_linpos = sub2ind(...
            size( gradc_unsc ), bimagVarCoilIdxShift, constCoilIdxShift  );
        bimag_minmax_sign_vec = repmat(...
            sign( bimag_minmax_arr( bimag_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_bimag_linpos ) = ...
            ( bimag_minmax_sign_vec .* shapefunvals_bimag ) .* bimag_sc( bimagVarCoilLocalIdxShift );

    end    

    gradc = gradc_unsc / RFMaxConstr;

end

end