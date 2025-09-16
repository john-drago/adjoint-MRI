function [ waveform ] = generateChebPlotWaveform_base( p, opt )
% [ RF, Grad, Shim, RFphasor, Freq, tvec, params ] fields of waveform
% This function will generate a waveform that can be plotted based on the
% unscaled variables in pSc.

% Get unscaled parameters -- currently expecting unscaled parameters to be
% passed
% scVec = opt.scVec;
% p = pSc( : ) .* scVec;

% Initialize Waveforms
numZCoils = opt.numZCoils;
numXYCoils = opt.numXYCoils;

tdom = opt.tdom;

% Generate time vector
pulseLength = opt.pulseLength;
sigDigitTime = 6;

dtInterp = opt.dt / 4;
tvecInterp = round( ( 0 : dtInterp : pulseLength ), sigDigitTime, 'significant' );
dtvecInterp = dtInterp * ones( size( tvecInterp ) );
dtvecInterp( 1 ) = 0;
dtvecInterp( end ) = 0;

tvecOpt = opt.opt_tvec;
tvecOptPad = [ 0, tvecOpt, pulseLength ];

numTimePoints = length( tvecInterp );
% dtError = eps( 1e3 );

% RF = complex( zeros( numXYCoils, ( numTimePoints ) ) );
% Grad = zeros( 3, ( numTimePoints ) );

%% Deal with Freq
freqInterpMethod = 'nearest';
dwxyvec_optpad = [ opt.dwxyvec(1), opt.dwxyvec, opt.dwxyvec(end) ];
Freq = interpolateWaveforms( tvecOptPad, dwxyvec_optpad / ( 2*pi ), tvecInterp, freqInterpMethod );
dwxyvec = ( 2 * pi ) * Freq;

%% Add RF points to RF array
numCheb_RF = opt.numCheb_RF;

% breal
breal_opt_idx = opt.breal_idx;
breal_opt_vec = p( breal_opt_idx );
breal_rshp = reshape( breal_opt_vec, [ numCheb_RF, numXYCoils ] );

% bimag
bimag_opt_idx = opt.bimag_idx;
bimag_opt_vec = p( bimag_opt_idx );
bimag_rshp = reshape( bimag_opt_vec, [ numCheb_RF, numXYCoils ] );

% bcomp
breal_wv = transpose( evalChebClenshaw( tvecInterp, breal_rshp, tdom ) );
bimag_wv = transpose( evalChebClenshaw( tvecInterp, bimag_rshp, tdom ) );

bcomp_wv = breal_wv + 1j * bimag_wv;

RF = bcomp_wv;
RFphasor = RF;

%% Add Grad points to Grad array
numCheb_grad = opt.numCheb_grad;

% grad
grad_opt_idx = opt.grad_idx;
grad_opt_vec = p( grad_opt_idx );
grad_rshp = reshape( grad_opt_vec, [ numCheb_grad, 3 ] );

Grad = transpose( evalChebClenshaw( tvecInterp, grad_rshp, tdom ) );

%% Add shim points to shim array
if numZCoils > 0
    
    numCheb_shim = opt.numCheb_shim;

    shim_opt_idx = opt.shim_idx;
    shim_opt_vec = p( shim_opt_idx );
    shim_rshp = reshape( shim_opt_vec, [ numCheb_shim, numZCoils ] );
    
    Shim = transpose( evalChebClenshaw( tvecInterp, shim_rshp, tdom ) );

else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_rshp = breal_rshp;
waveform.bimag_rshp = bimag_rshp;
waveform.numCheb_RF = numCheb_RF;
waveform.opt_tvec = tvecOpt;

waveform.Grad = Grad;
waveform.grad_rshp = grad_rshp;
waveform.numCheb_grad = numCheb_grad;

waveform.Shim = Shim;
if numZCoils > 0
    waveform.shim_rshp = shim_rshp;
    waveform.numCheb_shim = numCheb_shim;
end

waveform.numZCoils = numZCoils;
waveform.numXYCoils = numXYCoils;

waveform.numPos = opt.numPos;
waveform.numVars = opt.numVars;

waveform.dwxyvec = dwxyvec;
waveform.constantRotatingFrame = opt.constantRotatingFrame;

waveform.numTimePoints = numTimePoints;
waveform.tvec = tvecInterp;
waveform.dt = dtInterp;
waveform.dtvec = dtvecInterp;
waveform.tdom = tdom;
waveform.pulseLength = pulseLength;

% waveform.Tn = opt.Tn;

waveform.Z0 = opt.Z0;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %