function startStop = processStartStopPPTiming( timing )

if isempty( timing )
    startStop = zeros( 0, 2 );
else
    % determine size of the timing array
    if ( size( timing, 1 ) == 1 ) || ( size( timing, 2 ) == 1 )
        timing = timing( : );

        Ntiming = length( timing );

        startStop = zeros( Ntiming- 1, 2 );

        for nn = 1:( Ntiming-1 )
            startStop( nn, : ) = [ timing( nn ), timing( nn+1 ) ];
        end

    elseif ( size( timing, 1 ) == 2 ) || ( size( timing, 2 ) == 2 )
        if ( size( timing, 1 ) == 2 )
            timing = transpose( timing );
        end
        startStop = timing;
    end

end

numDigAfterDecimal = 8;
startStop = round( startStop, numDigAfterDecimal );

end