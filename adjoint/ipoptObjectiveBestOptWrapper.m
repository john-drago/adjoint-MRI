function [ cost, gradcost, wrapoutput ] = ipoptObjectiveBestOptWrapper(...
    pSc, objectiveFn, funcs, options )

persistent iter bestIter bestCost pScOpt;

constrTolSave = 5.0e-2;

if ( nargout >=0 ) && ( nargout < 3 )
    
    pSc = pSc( : );
    
    if nargout < 2
        cost = objectiveFn( pSc );
    else
        [ cost, gradcost ] = objectiveFn( pSc );
    end

    if isempty( iter )
        iter = uint32( 0 );
    end
    iter = iter + 1;

    if isempty( bestCost )
        bestCost = inf;
    else

        if cost < bestCost
            boundsBin = all( ( pSc >= ( options.lb - constrTolSave ) ) &...
                ( pSc <= ( options.ub + constrTolSave ) ) );
            if boundsBin
                c = funcs.constraints( pSc );
                constBin = all( ( c >= ( options.cl - constrTolSave ) ) &...
                    ( c <= ( options.cu + constrTolSave ) ) );
                if constBin
                    pScOpt = pSc;
                    bestIter = iter;
                    bestCost = cost;
                end
            end
        end

    end

elseif nargout > 2

    wrapoutput = struct;

    cost = [];
    gradcost = [];
    wrapoutput.numIter = iter;
    wrapoutput.bestIter = bestIter;
    wrapoutput.bestCost = bestCost;
    wrapoutput.pScOpt = pScOpt;
end

end