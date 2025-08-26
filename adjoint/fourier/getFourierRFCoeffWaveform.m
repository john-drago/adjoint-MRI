function c = getFourierRFCoeffWaveform( tvec, RFwvform, orderRF )

if all( RFwvform == 0 )
    c = 0;
    return
end

%% Assume waveform is in standard first dimension is waveform and second dimension is timepoints
RFwvform = transpose( RFwvform );
tvec = tvec( : );
tdom = [ tvec( 1 ), tvec( end ) ];

%% Loop over each waveform
orderRF = double( orderRF );
ncoeff = min( ( orderRF * 10 ), length( tvec ) );

wvformsampfun = makeWaveformSampleFunction( tvec, RFwvform );
[ wvformcoeffs, ~ ] = fourierCoeffFromFun( wvformsampfun, tdom, ncoeff );
% [ wvformcoeffs, fveccoeffs ] = fourierCoeffFromFun( wvformsampfun, tdom, ncoeff );

zeroIdx = floor( ncoeff / 2 ) + 1;

c = wvformcoeffs( ( ( zeroIdx - orderRF ) : 1 : ( zeroIdx + orderRF ) ), : );

end

