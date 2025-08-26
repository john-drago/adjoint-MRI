function c = getFourierCoeffWaveform( tvec, wvform, orderMax )

if all( wvform == 0 )
    c = 0;
    return
end

%% Assume waveform is in standard first dimension is waveform and second dimension is timepoints
wvform = transpose( wvform );
tvec = tvec( : );
tdom = [ tvec( 1 ), tvec( end ) ];

%% Loop over each waveform
orderMax = double( orderMax );
ncoeff = min( ( orderMax * 10 ), length( tvec ) );

wvformsampfun = makeWaveformSampleFunction( tvec, wvform );
[ wvformcoeffs, ~ ] = fourierCoeffFromFun( wvformsampfun, tdom, ncoeff );
% [ wvformcoeffs, fveccoeffs ] = fourierCoeffFromFun( wvformsampfun, tdom, ncoeff );

zeroIdx = floor( ncoeff / 2 ) + 1;
widthIdx = floor( ncoeff / 2 ) - 1;

posIdx = ( zeroIdx + 1 ) : 1 : ( zeroIdx + widthIdx );
negIdx = ( zeroIdx - 1 ) : -1 : ( zeroIdx - widthIdx );

coscoeff = real( [...
    wvformcoeffs( zeroIdx, : );...
    wvformcoeffs( posIdx, : ) + wvformcoeffs( negIdx, : );...
    ] );

sincoeff = real( [...
    wvformcoeffs( zeroIdx, : );...
    1j * ( wvformcoeffs( posIdx, : ) - wvformcoeffs( negIdx, : ) );...
    ] );

c = [...
    coscoeff( 1:(orderMax+1), : );...
    sincoeff( 2:(orderMax+1), : );...
    ];

end

