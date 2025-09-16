function [ peakEH10m, peakAvgEH10m, EH10m ] = calcEH10m( ...
    VOPs, RFphasor, dtvec, pulseLength, dutyCycle )

numVOPs = size( VOPs, 3 );
numTimePoints = size( dtvec, 2 );

EH10m = zeros( size( VOPs, 3 ), size( dtvec, 2 ) );

for tt = 1:numTimePoints
    RFconjVOP = repmat( ctranspose( RFphasor( :, tt ) ), [ 1, 1, numVOPs ] );
    RFVOP = repmat( RFphasor( :, tt ), [ 1, 1, numVOPs ] );

    EH10m( :, tt ) = sqrt( real( squeeze( pagemtimes( RFconjVOP, pagemtimes( VOPs, RFVOP ) ) ) ) );
end

peakEH10m = max( EH10m, [], 'all' );
avgEH10m = ( ( dutyCycle ) / ( pulseLength ) ) * ( EH10m * transpose( dtvec ) ); 
peakAvgEH10m = max( avgEH10m );
end