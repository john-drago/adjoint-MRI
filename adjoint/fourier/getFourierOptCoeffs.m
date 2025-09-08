function c = getFourierOptCoeffs( wv, orderRF, orderGrad, orderShim )

if ~isfield( wv, 'RF' )
    wv.RF = complex( wv.breal, wv.bimag );
end

%% Process RF waveforms
if ~isempty( orderRF )
    c.bcomp_coeffs = getFourierRFCoeffWaveform( wv.tvec, wv.RF, orderRF );
    c.breal_coeffs = real( c.bcomp_coeffs );
    c.bimag_coeffs = imag( c.bcomp_coeffs );
else
    c.bcomp_coeffs = [];
    c.breal_coeffs = [];
    c.bimag_coeffs = [];
end

%% Process Grad waveforms
if ~isempty( orderGrad )
    c.grad_coeffs = getFourierCoeffWaveform( wv.tvec, wv.Grad, orderGrad );
else
    c.grad_coeffs = [];
end

%% Process Shim waveforms
if ~isempty( orderShim )
    c.shim_coeffs = getChebCoeffWaveform( wv.tvec, wv.Shim, orderShim );
else
    c.shim_coeffs = [];
end

end