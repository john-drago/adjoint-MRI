function [ totalRFPower, maxRFPower, pulsePower_timepoint ] = calcTotalMaxRFPower( ...
    RFphasor, dtvec, pulseLength, Z0, dutyCycle )

pulsePower_timepoint = ( 1 / ( 2 * Z0 ) ) * abs( RFphasor ).^2;
pulsePower_avg = ( dutyCycle / pulseLength ) * sum( dtvec .* pulsePower_timepoint, 2 );

totalRFPower = sum( pulsePower_avg, 1 );
maxRFPower = max( pulsePower_avg, [], 'all' );

end