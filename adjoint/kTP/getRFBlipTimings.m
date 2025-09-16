function [ RFTiming, RFCenters, blipTiming, blipCenters ] = ...
    getRFBlipTimings( num_kTP, RFLength, blipLength )

timeResDigits = 6; % round to microseconds
% Assign timings
RFTiming = zeros( num_kTP, 2 );
RFCenters = zeros( num_kTP, 1 );
blipTiming = zeros( (num_kTP-1), 2 );
blipCenters = zeros( (num_kTP-1), 1 );

for nn = 1:num_kTP
    RFTiming( nn, 1 ) = ( nn - 1 ) * ( RFLength + blipLength);
    RFTiming( nn, 2 ) = RFLength + ( nn - 1 ) * ( RFLength + blipLength);

    RFCenters( nn, 1 ) = 0.5 * RFLength + ( nn - 1 ) * ( RFLength + blipLength);
    
    if nn ~= num_kTP
        blipTiming( nn, 1 ) = nn * RFLength + ( nn - 1 ) * blipLength;
        blipTiming( nn, 2 ) = nn * RFLength + ( nn ) * blipLength;

        blipCenters( nn, 1 ) = nn * RFLength + ( nn - 0.5 ) * blipLength;
    end
end

RFTiming = round( RFTiming, timeResDigits );
blipTiming = round( blipTiming, timeResDigits );

end