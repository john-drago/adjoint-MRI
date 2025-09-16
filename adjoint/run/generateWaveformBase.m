function wv = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt)

wv = struct;

wv.gyro = opt.gyro;

wv.breal = breal;
wv.bimag = bimag;
wv.brealphasor = brealphasor;
wv.bimagphasor = bimagphasor;

wv.Grad = Grad;

wv.Shim = Shim;

wv.numZCoils = opt.numZCoils;
wv.numXYCoils = opt.numXYCoils;

wv.numPos = opt.numPos;
wv.numVars = opt.numVars;

wv.dutyCycle = opt.dutyCycle;

wv.numTimePoints = opt.numTimePoints;
wv.tvec = tvec;
wv.dtvec = dtvec;
wv.dwxyvec = dwxyvec;

wv.Z0 = opt.Z0;

wv.pulseLength = ( tvec( end ) + dtvec( end ) ) - ( tvec( 1 ) - dtvec( 1 ) );

if ~isfield( opt, 'constantRotatingFrame' )
    wv.constantRotatingFrame = true;
else
    wv.constantRotatingFrame = opt.constantRotatingFrame;
end

end