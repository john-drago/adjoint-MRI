function [ stop ] = checkMaxTimeFMin( ~, ~, state, maxTime )

stop = false;

persistent iter optTimeStart

switch state
    case 'init'
        iter = uint32( 0 );
        optTimeStart = datetime;

    case 'iter'
        iter = iter + 1;

        currOptDuration = seconds( datetime - optTimeStart );

        if currOptDuration >= maxTime
            stop = true;
            fprintf( "\nOptimization Terminated Due to Max Time Constraints\nFinal Opt Time:\t%.3f seconds\n", currOptDuration );
        end
        
end


end