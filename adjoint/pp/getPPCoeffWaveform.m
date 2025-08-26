function c = getPPCoeffWaveform( tvec, wvform, projSt, tol )

if nargin < 4
    tol = 1e-10;
end

if all( wvform == 0 )
    c = 0;
    return
end

%% Assume waveform is in standard first dimension is waveform and second dimension is timepoints
wvform = transpose( wvform );
tvec = tvec( : );

numwv = size( wvform, 2 );

%% Initialize coefficient list
c = zeros( projSt.numVarsPerChannel, numwv );

for ww = 1:numwv
    c( :, ww ) = ppSFFromFun( wvform( :, ww ), tvec, projSt, tol );
end


end

