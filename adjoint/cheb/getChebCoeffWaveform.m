function c = getChebCoeffWaveform( tvec, wvform, tol )

if all( wvform == 0 )
    c = 0;
    return
end

if nargin < 3
    tol = 1e-3;
end

%% Assume waveform is in standard first dimension is waveform and second dimension is timepoints
wvform = transpose( wvform );
tvec = tvec( : );
tdom = [ tvec( 1 ), tvec( end ) ];

%% Loop over each waveform
nw = size( wvform, 2 );
wvformcoeffs = cell( 1, nw );
wvformsize = zeros( 1, nw );

for nn = 1:nw
    wvformsampfun = makeWaveformSampleFunction( tvec, wvform( :, nn ) );
    wvformcoeffs{ 1, nn } = chebCoeffFromFun( wvformsampfun, tdom, tol );
    wvformsize( nn ) = size( wvformcoeffs{ 1, nn }, 1 );
end

ncmax = max( wvformsize );

%% Make all coefficient vectors the same size
for nn = 1:nw
    wvformcoeffs{ 1, nn } = [ wvformcoeffs{ 1, nn }; zeros( ( ncmax - wvformsize( nn ) ), 1 ) ];
end

%% Return coefficient vector
cfull = cell2mat( wvformcoeffs );
absc = abs( cfull );

%% Try to trim coefficients
cenergy = sum( absc, 1 );
cenergyNonZero = cenergy > 0;
ccumsum = cumsum( absc( :, cenergyNonZero )./cenergy( :, cenergyNonZero ), 1 );
cnumtol = sum( ccumsum < ( 1 - tol ), 1 );

c = cfull( 1:min( cnumtol ), : );

end

