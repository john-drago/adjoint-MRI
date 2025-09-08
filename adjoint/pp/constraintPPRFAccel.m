function [ c, gradc ] = constraintPPRFAccel( pSc, opt )

RFAccelConstr = opt.RFAccel_constr;

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
brealcheb_arr = permute( reshape( brealcheb_rshp,...
    [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] ), [ 1, 3, 2 ] );
ddbrealcheb_arr = pagemtimes( opt.D2_RF, brealcheb_arr );

ddbreal_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
ddbreal_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
ddbreal_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
ddbreal_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numVarsPerChannel_RF, numXYCoils ] );
bimagcheb_rshp = varsToChebByPeriods_RF * bimag_rshp;
bimagcheb_arr = permute( reshape( bimagcheb_rshp,...
    [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] ), [ 1, 3, 2 ] );
ddbimagcheb_arr = pagemtimes( opt.D2_RF, bimagcheb_arr );

ddbimag_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
ddbimag_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
ddbimag_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
ddbimag_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

for pp = 1:numPPPeriods_RF
    [ ddbreal_minmax_arr( :, :, pp ), ddbreal_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( ddbrealcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );
    [ ddbimag_minmax_arr( :, :, pp ), ddbimag_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( ddbimagcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );

    [ ddbreal_absmax_arr( :, pp ), ddbreal_absmax_rowpos_arr( :, pp ) ] = max( abs( ddbreal_minmax_arr( :, :, pp ) ), [], 1 );
    [ ddbimag_absmax_arr( :, pp ), ddbimag_absmax_rowpos_arr( :, pp ) ] = max( abs( ddbimag_minmax_arr( :, :, pp ) ), [], 1 );
end

ddbreal_minmax_arr = permute( ddbreal_minmax_arr, [ 1, 3, 2 ] );
ddbimag_minmax_arr = permute( ddbimag_minmax_arr, [ 1, 3, 2 ] );

ddbreal_minmax_pos_arr = permute( ddbreal_minmax_pos_arr, [ 1, 3, 2 ] );
ddbimag_minmax_pos_arr = permute( ddbimag_minmax_pos_arr, [ 1, 3, 2 ] );

ddbreal_absmax_rowpos_arr = transpose( ddbreal_absmax_rowpos_arr );
ddbimag_absmax_rowpos_arr = transpose( ddbimag_absmax_rowpos_arr );

ddbreal_absmax_arr = transpose( ddbreal_absmax_arr );
ddbimag_absmax_arr = transpose( ddbimag_absmax_arr );

ddbreal_absmax = reshape( ddbreal_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );
ddbimag_absmax = reshape( ddbimag_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );

c_unsc = [ ddbreal_absmax; ddbimag_absmax ] - RFAccelConstr;
c = c_unsc / RFAccelConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2 * numXYCoils * numPPPeriods_RF );

    numXYCoils_numPPPeriods_arr = repmat(  1:numXYCoils, [ numPPPeriods_RF, 1 ] );
    numPPPeriods_numXYCoils_arr = repmat( transpose( 1:numPPPeriods_RF ), [ 1, numXYCoils ] );

    % breal
    ddbreal_absmax_linpos = sub2ind( size( ddbreal_minmax_arr ),...
        ddbreal_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    ddbreal_absmax_pos = ddbreal_minmax_pos_arr( ddbreal_absmax_linpos );

    % bimag
    ddbimag_absmax_linpos = sub2ind( size(ddbimag_minmax_arr ),...
        ddbimag_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    ddbimag_absmax_pos = ddbimag_minmax_pos_arr( ddbimag_absmax_linpos );

    for pp = 1:numPPPeriods_RF

        % breal
        ddbreal_minmax_Tn_samples = ...
            evalChebClenshaw( ddbreal_absmax_pos( pp, : ), eye( numPPShape_RF-2 ), PPStartStop_RF( pp, : ) );

        bVarIdxsLength = opt.PPVarIdxs_RF(pp,2)-opt.PPVarIdxs_RF(pp,1)+1;
        brealVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        brealVarCoilIdxShift = opt.breal_idx( brealVarCoilLocalIdxShift );
        constCoilIdxShift = repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );

        shapefunvals_breal = transpose( ddbreal_minmax_Tn_samples * ( opt.D2_RF( :, :, pp ) * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) ) ) );

        gradc_breal_linpos = sub2ind(...
            size( gradc_unsc ), brealVarCoilIdxShift, constCoilIdxShift  );
        breal_minmax_sign_vec = repmat(...
            sign( ddbreal_minmax_arr( ddbreal_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_breal_linpos ) = ...
            ( breal_minmax_sign_vec .* shapefunvals_breal ) .* breal_sc( brealVarCoilLocalIdxShift );

        % bimag
        bimag_minmax_Tn_samples = ...
            evalChebClenshaw( ddbimag_absmax_pos( pp, : ), eye( numPPShape_RF-2 ), PPStartStop_RF( pp, : ) );

        bimagVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        bimagVarCoilIdxShift = opt.bimag_idx( bimagVarCoilLocalIdxShift );
        constCoilIdxShift = numPPPeriods_RF * numXYCoils +...
            repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );

        shapefunvals_bimag = transpose( bimag_minmax_Tn_samples * ( opt.D2_RF( :, :, pp ) * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) ) ) );

        gradc_bimag_linpos = sub2ind(...
            size( gradc_unsc ), bimagVarCoilIdxShift, constCoilIdxShift  );
        bimag_minmax_sign_vec = repmat(...
            sign( ddbimag_minmax_arr( ddbimag_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_bimag_linpos ) = ...
            ( bimag_minmax_sign_vec .* shapefunvals_bimag ) .* bimag_sc( bimagVarCoilLocalIdxShift );

    end

    gradc = gradc_unsc / RFAccelConstr;

end

end