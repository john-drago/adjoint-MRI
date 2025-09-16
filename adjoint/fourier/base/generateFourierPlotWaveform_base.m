function [ waveform ] = generateFourierPlotWaveform_base( p, opt )
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
tvecOptPad = [ 0, tvecOpt, pulseLength ];

numTimePoints = length( tvecInterp );

omega0 = opt.omega0;

%% Deal with Freq
freqInterpMethod = 'nearest';
dwxyvec_optpad = [ opt.dwxyvec(1), opt.dwxyvec, opt.dwxyvec(end) ];
Freq = interpolateWaveforms( tvecOptPad, dwxyvec_optpad / ( 2*pi ), tvecInterp, freqInterpMethod );
dwxyvec = ( 2 * pi ) * Freq;

%% Add RF points to RF array
numFourier_RF = opt.numFourier_RF;

% b opt
breal_opt_idx = opt.breal_idx;
breal_opt_vec = p( breal_opt_idx );
breal_opt_rshp = reshape( breal_opt_vec, [ numFourier_RF, numXYCoils ] );

bimag_opt_idx = opt.bimag_idx;
bimag_opt_vec = p( bimag_opt_idx );
bimag_opt_rshp = reshape( bimag_opt_vec, [ numFourier_RF, numXYCoils ] );

bcomp_opt_rshp = complex(...
    breal_opt_rshp, bimag_opt_rshp );

% bcomp
FBRF = getFourierRFBasisMatrix( opt.orderFourier_RF, tvecInterp, omega0 );

bcomp_wv = FBRF * bcomp_opt_rshp;

bcomp_wv = transpose( bcomp_wv );

RF = bcomp_wv;
RFphasor = RF;

%% Add Grad points to Grad array
numFourier_grad = opt.numFourier_grad;

% grad
grad_opt_idx = opt.grad_idx;
grad_opt_vec = p( grad_opt_idx );
grad_opt_rshp = reshape( grad_opt_vec, [ numFourier_grad, 3 ] );

FBgrad = getFourierBasisMatrix( opt.orderFourier_grad, tvecInterp, omega0 );

Grad = FBgrad * grad_opt_rshp;

Grad = transpose( Grad );

%% Add shim points to shim array
if numZCoils > 0
    
    numFourier_shim = opt.numFourier_shim;

    % shim
    shim_opt_idx = opt.shim_idx;
    shim_opt_vec = p( shim_opt_idx );
    shim_opt_rshp = reshape( shim_opt_vec, [ numFourier_shim, numZCoils ] );

    FBshim = getFourierBasisMatrix( opt.orderFourier_shim, tvecInterp, omega0 );

    Shim = FBshim * shim_opt_rshp;

    Shim = transpose( Shim );

else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_rshp = breal_opt_rshp;
waveform.bimag_rshp = bimag_opt_rshp;
waveform.bcomp_rshp = bcomp_opt_rshp;
waveform.numFourier_RF = numFourier_RF;
waveform.opt_tvec = tvecOpt;

waveform.Grad = Grad;
waveform.grad_rshp = grad_opt_rshp;
waveform.numFourier_grad = numFourier_grad;

waveform.Shim = Shim;
if numZCoils > 0
    waveform.shim_rshp = shim_opt_rshp;
    waveform.numFourier_shim = numFourier_shim;
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

% waveform.FBRF = opt.FBRF;
% waveform.FBgrad = opt.FBgrad;
% waveform.FBshim = opt.FBshim;

waveform.Z0 = opt.Z0;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %