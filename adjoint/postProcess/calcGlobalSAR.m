function [ peakGlobalSAR, avgGlobalSAR, globalSAR ] = calcGlobalSAR( ...
    QGlobal, RFphasor, dtvec, pulseLength, dutyCycle )

numTimePoints = size( dtvec, 2 );

globalSAR = zeros( 1, size( dtvec, 2 ) );

for tt = 1:numTimePoints
    globalSAR( 1, tt ) = real( RFphasor( :, tt )' * QGlobal * RFphasor( :, tt ) );   
end

peakGlobalSAR = max( globalSAR );
avgGlobalSAR = ( ( dutyCycle ) / ( pulseLength ) ) * ( globalSAR * dtvec.' );

end