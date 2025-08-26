function wvsampfun = makeWaveformSampleFunction( tvec, wvform, method, extrap )
tvec = tvec( : );

if size( wvform, 1 ) == length( tvec )
elseif size( wvform, 2 ) == length( tvec )
    wvform = transpose( wvform );
else
    error( "Size of waveform not compatible for interpolation." )
end

if nargin < 3
    method = 'pchip';
end
if nargin < 4
    extrap = nan;
end

wvsampfun = @( tsamp ) interp1( tvec, wvform, tsamp( : ), method, extrap );

end