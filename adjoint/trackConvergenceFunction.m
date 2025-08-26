function [ stop, fval_final, optTime_final, funccount_final ] = trackConvergenceFunction(...
    ~, optimValues, state )

stop = false;

persistent iter optTimeStart fval_iter optTime_iter funccount_iter;

switch state
    case 'init'
        iter = uint32( 0 );

    case 'iter'
        iter = iter + 1;

        if iter == 1
            optTimeStart = datetime;
            optTime_iter( iter, 1 ) = 0;
        else
            optTime_iter( iter, 1 ) = seconds( datetime - optTimeStart );
        end
        fval_iter( iter, 1 ) = optimValues.fval;
        funccount_iter( iter, 1 ) = optimValues.funccount;

    case 'done'

        optTime_iter = optTime_iter( : );
        fval_iter = fval_iter( : );
        funccount_iter = funccount_iter( : );

    case 'out'

        optTime_final = optTime_iter( 1:iter, 1 );
        fval_final = fval_iter( 1:iter, 1 );
        funccount_final = funccount_iter( 1:iter, 1 );

end

end