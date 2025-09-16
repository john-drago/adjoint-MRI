function [ cost, gradcost, outStruct ] =...
    costFnTrackDecisionVariablesWrapper( pSc, costFn, const, spacingTrackDecisionVariables, constrTolSave )

if nargin < 5
    constrTolSave = 5.0e-2;
end

persistent iter idxCtr bestCost optTimeStart...
    fval_iters optTime_iters funccount_iters...
    fvalraw_iters const_iters...
    xtrack_iters bestpSc bestIter;

if ( nargout >=0 ) && ( nargout < 3 )

    if nargout > 1
        [ cost, gradcost ] = costFn( pSc );
    else
        [ cost ] = costFn( pSc );
    end

    if isempty( iter )
        iter = uint32( 0 );
    end
    iter = iter + uint32( 1 );
    
    if isempty( optTimeStart )
        optTimeStart = datetime;
    end

    if isempty( bestCost )
        bestCost = inf;
    end

    [ constBin, maxConstrViol ] = ( runConstraintsConvergence( pSc, const, constrTolSave ) );

    if ( cost < bestCost )
        if constBin || ( iter == uint32( 1 ) )
            saveBest = true;
            bestCost = cost;
            bestIter = iter;
            bestpSc = pSc( : );
        else
            saveBest = false;
        end
    else
        saveBest = false;
    end
    
    spacingBin = ( ( uint32( mod( iter, spacingTrackDecisionVariables ) ) == uint32( 0 ) ) || ...
            ( iter == uint32( 1 ) ) );

    if ( saveBest || spacingBin )

        if isempty( idxCtr )
            idxCtr = uint32( 0 );
        end
        idxCtr = idxCtr + uint32( 1 );

        if isempty( xtrack_iters )
            xtrack_iters = [];
        end

        if saveBest
            fval_iters( idxCtr, 1 ) = bestCost;
            xtrack_iters( :, idxCtr ) = bestpSc;
        elseif constBin
            fval_iters( idxCtr, 1 ) = cost;
            xtrack_iters( :, idxCtr ) = pSc( : );
        else
            fval_iters( idxCtr, 1 ) = bestCost;
            xtrack_iters( :, idxCtr ) = bestpSc;
        end

        fvalraw_iters( idxCtr, 1 ) = cost;
        const_iters( idxCtr, 1 ) = maxConstrViol;

        funccount_iters( idxCtr, 1 ) = iter;

        if idxCtr > uint32( 1 )
            optTime_iters( idxCtr, 1 ) = seconds( datetime - optTimeStart );
        else
            optTime_iters( idxCtr, 1 ) = 0;
        end

    end

elseif nargout > 2

    [ cost, gradcost ] = costFn( pSc );

    [ constBin, maxConstrViol ] = ( runConstraintsConvergence( pSc, const, constrTolSave ) );

    if ( cost < bestCost ) && ( constBin )

        iter = iter + uint32( 1 );
        idxCtr = idxCtr + uint32( 1 );

        bestCost = cost;
        bestIter = iter;
        bestpSc = pSc( : );
        
        xtrack_iters( :, idxCtr ) = bestpSc;
        fval_iters( idxCtr, 1 ) = bestCost;
        optTime_iters( idxCtr, 1 ) = seconds( datetime - optTimeStart );
        funccount_iters( idxCtr, 1 ) = iter;
        fvalraw_iters( idxCtr, 1 ) = cost;
        const_iters( idxCtr, 1 ) = maxConstrViol;

    end

    outStruct.fval_vars_iters = fval_iters( : );
    outStruct.optTime_vars_iters = optTime_iters( : );
    outStruct.funccount_vars_iters = funccount_iters( : );
    outStruct.fvalraw_vars_iters = fvalraw_iters( : );
    outStruct.const_vars_iters = const_iters( : );
    outStruct.xtrack_vars_iters = xtrack_iters;

    outStruct.pScOpt = bestpSc( : );

    outStruct.bestIter = bestIter;
    outStruct.bestCost = bestCost;

end

end
