function c = getChebOptCoeffs( wv, tol )

if nargin < 2
    tol = 1e-4;
end

if ~isfield( wv, 'breal' ) || ~isfield( wv, 'bimag' )
    wv.breal = real( wv.RF );
    wv.bimag = imag( wv.RF );
end

%% Process RF waveforms
numRFwv = size( wv.breal, 1 );
c.breal_coeffs = getChebCoeffWaveform( wv.tvec, wv.breal, tol );
numbrealcoeff = size( c.breal_coeffs, 1 );
c.bimag_coeffs = getChebCoeffWaveform( wv.tvec, wv.bimag, tol );
numbimagcoeff = size( c.bimag_coeffs, 1 );

if numbrealcoeff <= numbimagcoeff
    c.breal_coeffs = [ c.breal_coeffs; zeros( numbimagcoeff-numbrealcoeff, numRFwv ) ];
else
    c.bimag_coeffs = [ c.bimag_coeffs; zeros( numbrealcoeff-numbimagcoeff, numRFwv ) ];
end

%% Process Grad waveforms
c.grad_coeffs = getChebCoeffWaveform( wv.tvec, wv.Grad, tol );

%% Process Shim waveforms
c.shim_coeffs = getChebCoeffWaveform( wv.tvec, wv.Shim, tol );

end