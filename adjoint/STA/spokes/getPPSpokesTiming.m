function pulse = getPPSpokesTiming( pulse, spokesTiming )

numSpokes = pulse.numSpokes;
centralSpokeLength = spokesTiming.centralSpokeLength;
nonCentralSpokeLength = spokesTiming.nonCentralSpokeLength;
rephaseLength = spokesTiming.rephaseLength;
lastBlipLength = spokesTiming.lastBlipLength;
spokeBlipLength = spokesTiming.spokeBlipLength;
pulseLength = spokesTiming.pulseLength;

rephasePeriod = [ pulseLength-rephaseLength, pulseLength ];
spokesPeriods = zeros( numSpokes, 2 );
spokesPeriods( end, : ) = [ rephasePeriod( 1 )-centralSpokeLength, rephasePeriod( 1 ) ];

if numSpokes > 1
    blipPeriods = zeros( numSpokes-1, 2 );
    for nn = (numSpokes-1):-1:1
        if nn == ( numSpokes - 1 )
            blipPeriods( nn, : ) = [ spokesPeriods(nn+1, 1)-lastBlipLength, spokesPeriods(nn+1,1) ];
        else
            blipPeriods( nn, : ) = [ spokesPeriods(nn+1, 1)-spokeBlipLength, spokesPeriods(nn+1,1) ];
        end
        spokesPeriods( nn, : ) = [ blipPeriods( nn, 1 )-nonCentralSpokeLength, blipPeriods( nn, 1 ) ];
    end
end

% set piecewise polynomial factors
if numSpokes > 1
    % pulse.timing_RF = zeros( 2 + ( numSpokes - 1 ), 2 );
    % pulse.timing_grad = zeros( 2 + ( numSpokes - 1 ), 2 );

    pulse.timing_RF = zeros( numSpokes + 1, 2 );
    pulse.timing_grad = zeros( 1 + numSpokes + ( numSpokes - 1 ), 2 );
else
    pulse.timing_RF = zeros( 2, 2 );
    pulse.timing_grad = zeros( 2, 2 );
end

for nn = numSpokes:-1:1

    if nn == numSpokes
        pulse.timing_RF( end, : ) = rephasePeriod( 1, : );
        pulse.timing_RF( (end-1), : ) = spokesPeriods( end, : );

        pulse.timing_grad( end, : ) = rephasePeriod;
        pulse.timing_grad( (end-1), : ) = spokesPeriods( end, : );

    % elseif nn == ( numSpokes - 1 )
    % 
    %     pulse.timing_RF( nn+1, : ) = blipPeriods( nn, : );
    %     pulse.timing_RF( nn, : ) = spokesPeriods( nn, : );
    % 
    %     pulse.timing_grad( nn+1, : ) = blipPeriods( nn, : );
    %     pulse.timing_grad( nn, : ) = spokesPeriods( nn, : );

    else
        % pulse.timing_RF( nn, : ) = [ spokesPeriods( nn, 1 ), blipPeriods( nn, 2 ) ];
        % pulse.timing_grad( nn, : ) = [ spokesPeriods( nn, 1 ), blipPeriods( nn, 2 ) ];

        pulse.timing_RF( nn, 1 ) = spokesPeriods( nn, 1 );
        pulse.timing_RF( nn, 2 ) = blipPeriods( nn, 2 );

        pulse.timing_grad( 2*nn-1, : ) = spokesPeriods( nn, : );
        pulse.timing_grad( 2*nn, : ) = blipPeriods( nn, : );
    end

end

numDigAfterDecimal = 10;
pulse.timing_RF = round( pulse.timing_RF, numDigAfterDecimal );
pulse.timing_grad = round( pulse.timing_grad, numDigAfterDecimal );

pulse.timing_shim = [];

end