function [ peakLocalSAR, peakAvgLocalSAR, localSAR ] = calcLocalSAR( ...
    VOPs, RFphasor, dtvec, pulseLength, dutyCycle )

numVOPs = size( VOPs, 3 );
numTimePoints = size( dtvec, 2 );

localSAR = zeros( size( VOPs, 3 ), size( dtvec, 2 ) );

for tt = 1:numTimePoints
    RFconjVOP = repmat( ctranspose( RFphasor( :, tt ) ), [ 1, 1, numVOPs ] );
    RFVOP = repmat( RFphasor( :, tt ), [ 1, 1, numVOPs ] );

    localSAR( :, tt ) = real( squeeze( pagemtimes( RFconjVOP, pagemtimes( VOPs, RFVOP ) ) ) );
end

peakLocalSAR = max( localSAR, [], 'all' );
avgLocalSAR = ( ( dutyCycle ) / ( pulseLength ) ) * ( localSAR * transpose( dtvec ) ); 
peakAvgLocalSAR = max( avgLocalSAR );
end