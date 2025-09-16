function [ c, gradc ] = constraintPPRFSlewRate( pSc, opt )

RFSlewRateConstr = opt.RFSlewRate_constr;

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
dbrealcheb_arr = pagemtimes( opt.D_RF, brealcheb_arr );

dbreal_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
dbreal_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
dbreal_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
dbreal_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numVarsPerChannel_RF, numXYCoils ] );
bimagcheb_rshp = varsToChebByPeriods_RF * bimag_rshp;
bimagcheb_arr = permute( reshape( bimagcheb_rshp,...
    [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] ), [ 1, 3, 2 ] );
dbimagcheb_arr = pagemtimes( opt.D_RF, bimagcheb_arr );

dbimag_minmax_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
dbimag_minmax_pos_arr = zeros( 2, numXYCoils, numPPPeriods_RF );
dbimag_absmax_arr = zeros( numXYCoils, numPPPeriods_RF );
dbimag_absmax_rowpos_arr = zeros( numXYCoils, numPPPeriods_RF );

for pp = 1:numPPPeriods_RF
    [ dbreal_minmax_arr( :, :, pp ), dbreal_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( dbrealcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );
    [ dbimag_minmax_arr( :, :, pp ), dbimag_minmax_pos_arr( :, :, pp ) ] = ...
        chebMinMax( dbimagcheb_arr( :, :, pp ), PPStartStop_RF( pp, : ) );

    [ dbreal_absmax_arr( :, pp ), dbreal_absmax_rowpos_arr( :, pp ) ] = max( abs( dbreal_minmax_arr( :, :, pp ) ), [], 1 );
    [ dbimag_absmax_arr( :, pp ), dbimag_absmax_rowpos_arr( :, pp ) ] = max( abs( dbimag_minmax_arr( :, :, pp ) ), [], 1 );
end

dbreal_minmax_arr = permute( dbreal_minmax_arr, [ 1, 3, 2 ] );
dbimag_minmax_arr = permute( dbimag_minmax_arr, [ 1, 3, 2 ] );

dbreal_minmax_pos_arr = permute( dbreal_minmax_pos_arr, [ 1, 3, 2 ] );
dbimag_minmax_pos_arr = permute( dbimag_minmax_pos_arr, [ 1, 3, 2 ] );

dbreal_absmax_rowpos_arr = transpose( dbreal_absmax_rowpos_arr );
dbimag_absmax_rowpos_arr = transpose( dbimag_absmax_rowpos_arr );

dbreal_absmax_arr = transpose( dbreal_absmax_arr );
dbimag_absmax_arr = transpose( dbimag_absmax_arr );

dbreal_absmax = reshape( dbreal_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );
dbimag_absmax = reshape( dbimag_absmax_arr, [ numPPPeriods_RF * numXYCoils, 1 ] );

c_unsc = [ dbreal_absmax; dbimag_absmax ] - RFSlewRateConstr;
c = c_unsc / RFSlewRateConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 2 * numXYCoils * numPPPeriods_RF );

    numXYCoils_numPPPeriods_arr = repmat(  1:numXYCoils, [ numPPPeriods_RF, 1 ] );
    numPPPeriods_numXYCoils_arr = repmat( transpose( 1:numPPPeriods_RF ), [ 1, numXYCoils ] );

    % breal
    dbreal_absmax_linpos = sub2ind( size( dbreal_minmax_arr ),...
        dbreal_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    dbreal_absmax_pos = dbreal_minmax_pos_arr( dbreal_absmax_linpos );

    % bimag
    dbimag_absmax_linpos = sub2ind( size(dbimag_minmax_arr ),...
        dbimag_absmax_rowpos_arr, numPPPeriods_numXYCoils_arr, numXYCoils_numPPPeriods_arr );
    dbimag_absmax_pos = dbimag_minmax_pos_arr( dbimag_absmax_linpos );

    for pp = 1:numPPPeriods_RF

        % breal
        dbreal_minmax_Tn_samples = ...
            evalChebClenshaw( dbreal_absmax_pos( pp, : ), eye( numPPShape_RF-1 ), PPStartStop_RF( pp, : ) );

        bVarIdxsLength = opt.PPVarIdxs_RF(pp,2)-opt.PPVarIdxs_RF(pp,1)+1;
        brealVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        brealVarCoilIdxShift = opt.breal_idx( brealVarCoilLocalIdxShift );
        constCoilIdxShift = repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );

        shapefunvals_breal = transpose( dbreal_minmax_Tn_samples * ( opt.D_RF( :, :, pp ) * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) ) ) );

        gradc_breal_linpos = sub2ind(...
            size( gradc_unsc ), brealVarCoilIdxShift, constCoilIdxShift  );
        breal_minmax_sign_vec = repmat(...
            sign( dbreal_minmax_arr( dbreal_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_breal_linpos ) = ...
            ( breal_minmax_sign_vec .* shapefunvals_breal ) .* breal_sc( brealVarCoilLocalIdxShift );

        % bimag
        bimag_minmax_Tn_samples = ...
            evalChebClenshaw( dbimag_absmax_pos( pp, : ), eye( numPPShape_RF-1 ), PPStartStop_RF( pp, : ) );

        bimagVarCoilLocalIdxShift = opt.PPVarIdxsCoilCell_RF{ pp };
        bimagVarCoilIdxShift = opt.bimag_idx( bimagVarCoilLocalIdxShift );
        constCoilIdxShift = numPPPeriods_RF * numXYCoils +...
            repmat( ( pp:numPPPeriods_RF:(numPPPeriods_RF*numXYCoils) ) , [ bVarIdxsLength, 1 ] );

        shapefunvals_bimag = transpose( bimag_minmax_Tn_samples * ( opt.D_RF( :, :, pp ) * opt.shapeFnChebCoeffs_RF( :, opt.PPVarSFIdxs_RF(pp,1):opt.PPVarSFIdxs_RF(pp,2) ) ) );

        gradc_bimag_linpos = sub2ind(...
            size( gradc_unsc ), bimagVarCoilIdxShift, constCoilIdxShift  );
        bimag_minmax_sign_vec = repmat(...
            sign( dbimag_minmax_arr( dbimag_absmax_linpos(pp,:) ) ), [ bVarIdxsLength, 1 ] );

        gradc_unsc( gradc_bimag_linpos ) = ...
            ( bimag_minmax_sign_vec .* shapefunvals_bimag ) .* bimag_sc( bimagVarCoilLocalIdxShift );

    end

    gradc = gradc_unsc / RFSlewRateConstr;

end

end