function [ waveform ] = generatePWCPlotWaveform_base( p, opt )
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

% Generate time vector
pulseLength = opt.pulseLength;
sigDigitTime = 6;

dtInterp = opt.dt / 4;
tvecInterp = round( ( 0 : dtInterp : pulseLength ), sigDigitTime, 'significant' );
dtvecInterp = dtInterp * ones( size( tvecInterp ) );
dtvecInterp( 1 ) = 0;
dtvecInterp( end ) = 0;

tvecOpt = opt.opt_tvec;
numTimePointsOpt = size( tvecOpt, 2 );
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
% breal
breal_opt_idx = opt.breal_idx;
breal_opt_vec = p( breal_opt_idx );
breal_opt = reshape( breal_opt_vec, [ numTimePointsOpt, numXYCoils ] ).';

% bimag
bimag_opt_idx = opt.bimag_idx;
bimag_opt_vec = p( bimag_opt_idx );
bimag_opt = reshape( bimag_opt_vec, [ numTimePointsOpt, numXYCoils ] ).';

% bcomp
bcomp_opt = breal_opt + 1j * bimag_opt;
bcomp_opt_pad = [ zeros(numXYCoils, 1), bcomp_opt, zeros(numXYCoils, 1) ];

RFphasor = interpolateWaveforms( tvecOptPad, bcomp_opt_pad, tvecInterp );
RF = RFphasor;

%% Add Grad points to Grad array
% grad
grad_opt_idx = opt.grad_idx;
grad_opt_vec = p( grad_opt_idx );
grad_opt = reshape( grad_opt_vec, [ numTimePointsOpt, 3 ] ).';
grad_opt_pad = [ zeros(3,1), grad_opt, zeros(3,1) ];

Grad = interpolateWaveforms( tvecOptPad, grad_opt_pad, tvecInterp );

%% Add shim points to shim array
if numZCoils > 0
    shim_opt_idx = opt.shim_idx;
    shim_opt_vec = p( shim_opt_idx );
    shim_opt = reshape( shim_opt_vec, [ numTimePointsOpt, numZCoils ] ).';
    shim_opt_pad = [ zeros(numZCoils,1), shim_opt, zeros(numZCoils,1) ];

    Shim = interpolateWaveforms( tvecOptPad, shim_opt_pad, tvecInterp );
else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_opt = breal_opt;
waveform.bimag_opt = bimag_opt;
waveform.bcomp_opt = bcomp_opt;
waveform.opt_tvec = tvecOpt;

waveform.Grad = Grad;
waveform.grad_opt = grad_opt;

waveform.Shim = Shim;
if numZCoils > 0
    waveform.shim_opt = shim_opt;
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
waveform.pulseLength = pulseLength;

waveform.Z0 = opt.Z0;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %