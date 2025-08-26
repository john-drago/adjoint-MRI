function [ c, gradc ] = constraintPPShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numZCoils = opt.numZCoils;
PPStartStop_shim = opt.PPStartStop_shim;
numVarsPerChannel_shim = opt.numVarsPerChannel_shim;
varsToChebByPeriods_shim = opt.varsToChebByPeriods_shim;
numPPShape_shim = opt.numPPShape_shim;
numPPPeriods_shim = opt.numPPPeriods_shim;
numDigitRound = 15;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );
shim_rshp = reshape( shim, [ numVarsPerChannel_shim, numZCoils ] );
shimcheb_rshp = varsToChebByPeriods_shim * shim_rshp;
shimcheb_arr = permute( reshape( shimcheb_rshp, [ numPPShape_shim, numPPPeriods_shim, numZCoils ] ), [ 1, 3, 2 ] );

currTotalMax_arr = zeros( numPPPeriods_shim, 1 );
currTotalPos_arr = zeros( numPPPeriods_shim, 1 );

absflag = true;

for pp = 1:numPPPeriods_shim
    shim_breakpts = chebRoots( shimcheb_arr( :, :, pp ), PPStartStop_shim( pp, : ) );
    shim_breakpts_sort = unique(...
        round( [ PPStartStop_shim( pp, 1 ); cell2mat( transpose( shim_breakpts ) ); PPStartStop_shim( pp, 2 ) ], numDigitRound ) );
    shim_breaks = chebRestrict( shimcheb_arr( :, :, pp ), shim_breakpts_sort, PPStartStop_shim( pp, : ), absflag );

    currTotalMax = 0;
    currTotalPos = 1;
    for bb = 1:size( shim_breaks, 1 )
        shim_bb = sum( cell2mat( shim_breaks( bb, : ) ), 2 );

        [ shim_bb_minmax, shim_bb_pos ] = chebMinMax( shim_bb, [ shim_breakpts_sort( bb ), shim_breakpts_sort( bb+1 ) ] );
        if shim_bb_minmax( 2 ) > currTotalMax
            currTotalMax = shim_bb_minmax( 2 );
            % currTotalIdx = bb;
            currTotalPos = shim_bb_pos( 2 );
        end
    end
    currTotalMax_arr( pp ) = currTotalMax;
    currTotalPos_arr( pp ) = currTotalPos;

end

c_unsc = currTotalMax_arr - shimTotalConstr;

c = c_unsc / shimTotalConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, numPPPeriods_shim );

    for pp = 1:numPPPeriods_shim

        shim_minmax_Tn_samples = ...
            evalChebClenshaw( currTotalPos_arr( pp ), eye( numPPShape_shim ), PPStartStop_shim( pp, : ) );

        shimVarIdxs = opt.PPVarIdxs_shim(pp,1):opt.PPVarIdxs_shim(pp,2);
        shimVarCoilLocalIdxShift = repmat( transpose( shimVarIdxs ), [ 1, numZCoils ] ) +...
            ((0:uint32(numZCoils-1))*(numVarsPerChannel_shim));
        shimVarCoilIdxShift = opt.shim_idx( shimVarCoilLocalIdxShift );
        
        shapefunvals_shim = transpose( shim_minmax_Tn_samples * ... 
            opt.shapeFnChebCoeffs_shim( :, opt.PPVarSFIdxs_shim(pp,1):opt.PPVarSFIdxs_shim(pp,2) ) );

        scshapefunvals_shim = transpose( shim_rshp( opt.PPVarIdxs_shim( pp, 1 ):opt.PPVarIdxs_shim( pp, 2 ), : ) )...
            * shapefunvals_shim;

        gradc_shim_linpos = sub2ind(...
            size( gradc_unsc ), shimVarCoilIdxShift, pp*ones( size( shimVarCoilIdxShift ) )  );
        gradc_unsc( gradc_shim_linpos ) = ...
            transpose( sign( scshapefunvals_shim ) ) .* ( shapefunvals_shim .* shim_sc( shimVarCoilLocalIdxShift ) );

    end

    gradc = gradc_unsc / shimTotalConstr;

end

end