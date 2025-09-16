function c = getPWCOptCoeffs( wv, tvecsamp )

if ~isfield( wv, 'breal' ) || ~isfield( wv, 'bimag' )
    wv.breal = real( wv.RF );
    wv.bimag = imag( wv.RF );
end

%% Process RF waveforms

brealsampfun = makeWaveformSampleFunction( wv.tvec, wv.breal );
c.breal = brealsampfun( tvecsamp( : ) );

bimagsampfun = makeWaveformSampleFunction( wv.tvec, wv.bimag );
c.bimag = bimagsampfun( tvecsamp( : ) );

%% Process Grad waveforms
gradsampfun = makeWaveformSampleFunction( wv.tvec, wv.Grad );
c.grad = gradsampfun( tvecsamp( : ) );

%% Process Shim waveforms
shimsampfun = makeWaveformSampleFunction( wv.tvec, wv.Shim );
c.shim = shimsampfun( tvecsamp( : ) );

end