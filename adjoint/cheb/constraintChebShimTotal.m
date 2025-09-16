function [ c, gradc ] = constraintChebShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

numZCoils = opt.numZCoils;
tdom = opt.tdom;
numDigitRound = 15;

shim_idx = opt.shim_idx;
shim_sc = opt.scVec( shim_idx );
shim = shim_sc .* pSc( shim_idx );

shim_sc_rshp = reshape( shim_sc, [ opt.numCheb_shim, numZCoils ] );
shim_idx_rshp = reshape( shim_idx, [ opt.numCheb_shim, numZCoils ] );
shim_rshp = reshape( shim, [ opt.numCheb_shim, numZCoils ] );

shim_breakpts = chebRoots( shim_rshp, tdom );
shim_breakpts_sort = unique( ...
     round( [ tdom( 1 ); cell2mat( transpose( shim_breakpts ) ); tdom( end ) ], numDigitRound ) );
absflag = true;
shim_breaks = chebRestrict( shim_rshp, shim_breakpts_sort, tdom, absflag );

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

c_unsc = currTotalMax - shimTotalConstr;

c = c_unsc / shimTotalConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, 1 );

    Tn = evalChebClenshaw( currTotalPos, eye( opt.numCheb_shim ), tdom );
    Tn_repmat = repmat( transpose( Tn ), [ 1, numZCoils ] );
    Tn_dot_shimcoeff = repmat( Tn * shim_rshp, [ opt.numCheb_shim, 1 ] );

    gradc_unsc( shim_idx_rshp ) =...
        sign( Tn_dot_shimcoeff ) .* Tn_repmat .* shim_sc_rshp;

    gradc = gradc_unsc / shimTotalConstr;

end

end