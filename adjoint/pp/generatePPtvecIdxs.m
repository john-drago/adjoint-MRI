function [ PPIdxs, numPPPeriods ] = generatePPtvecIdxs( PPStartStop, tvec, dttol )
if nargin < 3
    dttol = ( tvec( 2 ) - tvec( 1 ) ) / 10;
end

numPPPeriods = size( PPStartStop, 1 );
PPIdxs = zeros( numPPPeriods, 2, "uint32" );
for ii = 1:numPPPeriods
    PPIdxs( ii, 1 ) = find( tvec >= ( PPStartStop( ii, 1 ) - dttol  ), 1, 'first' );
    PPIdxs( ii, 2 ) = find( tvec <= ( PPStartStop( ii, 2 ) + dttol ), 1, 'last' );

    if ii > 1
        if PPIdxs( ii, 1 ) == PPIdxs( ii-1, 2 )
            PPIdxs( ii, 1 ) = PPIdxs( ii, 1 ) + 1;
        end
    end
end

end