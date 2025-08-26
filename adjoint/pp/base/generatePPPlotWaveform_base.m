function [ waveform ] = generatePPPlotWaveform_base( p, opt )
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
dtInterpTol = dtInterp / 10;

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
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;
numPPShape_RF = opt.numPPShape_RF;
numPPPeriods_RF = opt.numPPPeriods_RF;
PPStartStop_RF = opt.PPStartStop_RF;

% breal
breal_wv = zeros( numXYCoils, numTimePoints );
breal_opt_idx = opt.breal_idx;
breal_opt_vec = p( breal_opt_idx );
breal_rshp = reshape( breal_opt_vec, [ numVarsPerChannel_RF, numXYCoils ] );
breal_arr = reshape( opt.varsToChebByPeriods_RF * breal_rshp, [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] );
breal_arr = permute( breal_arr, [ 1, 3, 2 ] );

% bimag
bimag_wv = zeros( numXYCoils, numTimePoints );
bimag_opt_idx = opt.bimag_idx;
bimag_opt_vec = p( bimag_opt_idx );
bimag_rshp = reshape( bimag_opt_vec, [ numVarsPerChannel_RF, numXYCoils ] );
bimag_arr = reshape( opt.varsToChebByPeriods_RF * bimag_rshp, [ numPPShape_RF, numPPPeriods_RF, numXYCoils ] );
bimag_arr = permute( bimag_arr, [ 1, 3, 2 ] );

% bcomp
for pp = 1:numPPPeriods_RF
    timeIdx = ( ( tvecInterp ) >= ( PPStartStop_RF( pp, 1 ) - dtInterpTol ) ) &...
        ( ( tvecInterp ) <= ( PPStartStop_RF( pp, 2 ) + dtInterpTol ) );

    breal_wv( :, timeIdx ) = transpose( evalChebClenshaw(...
        tvecInterp( timeIdx ), breal_arr( :, :, pp ), PPStartStop_RF( pp, : ) ) );
    bimag_wv( :, timeIdx ) = transpose( evalChebClenshaw(...
        tvecInterp( timeIdx ), bimag_arr( :, :, pp ), PPStartStop_RF( pp, : ) ) );
end

bcomp_wv = complex( breal_wv, bimag_wv );

RF = bcomp_wv;
RFphasor = RF;

%% Add Grad points to Grad array
numVarsPerChannel_grad = opt.numVarsPerChannel_grad;
numPPShape_grad = opt.numPPShape_grad;
numPPPeriods_grad = opt.numPPPeriods_grad;
PPStartStop_grad = opt.PPStartStop_grad;

% grad
Grad = zeros( 3, numTimePoints );
grad_opt_idx = opt.grad_idx;
grad_opt_vec = p( grad_opt_idx );
grad_rshp = reshape( grad_opt_vec, [ numVarsPerChannel_grad, 3 ] );
grad_arr = reshape( opt.varsToChebByPeriods_grad * grad_rshp, [ numPPShape_grad, numPPPeriods_grad, 3 ] );
grad_arr = permute( grad_arr, [ 1, 3, 2 ] );

for pp = 1:numPPPeriods_grad
    timeIdx = ( ( tvecInterp ) >= ( PPStartStop_grad( pp, 1 ) - dtInterpTol ) ) &...
        ( ( tvecInterp ) <= ( PPStartStop_grad( pp, 2 ) + dtInterpTol ) );

    Grad( :, timeIdx ) = transpose( evalChebClenshaw(...
        tvecInterp( timeIdx ), grad_arr( :, :, pp ), PPStartStop_grad( pp, : ) ) );
    Grad( :, timeIdx ) = transpose( evalChebClenshaw(...
        tvecInterp( timeIdx ), grad_arr( :, :, pp ), PPStartStop_grad( pp, : ) ) );
end


%% Add shim points to shim array
if numZCoils > 0

    numVarsPerChannel_shim = opt.numVarsPerChannel_shim;
    numPPShape_shim = opt.numPPShape_shim;
    numPPPeriods_shim = opt.numPPPeriods_shim;
    PPStartStop_shim = opt.PPStartStop_shim;

    % shim
    Shim = zeros( numZCoils, numTimePoints );
    shim_opt_idx = opt.shim_idx;
    shim_opt_vec = p( shim_opt_idx );
    shim_rshp = reshape( shim_opt_vec, [ numVarsPerChannel_shim, numZCoils ] );
    shim_arr = reshape( opt.varsToChebByPeriods_shim * shim_rshp, [ numPPShape_shim, numPPPeriods_shim, numZCoils ] );
    shim_arr = permute( shim_arr, [ 1, 3, 2 ] );

    for pp = 1:numPPPeriods_shim
        timeIdx = ( ( tvecInterp ) >= ( PPStartStop_shim( pp, 1 ) - dtInterpTol ) ) &...
        ( ( tvecInterp ) <= ( PPStartStop_shim( pp, 2 ) + dtInterpTol ) );

        Shim( :, timeIdx ) = transpose( evalChebClenshaw(...
            tvecInterp( timeIdx ), shim_arr( :, :, pp ), PPStartStop_shim( pp, : ) ) );
        Shim( :, timeIdx ) = transpose( evalChebClenshaw(...
            tvecInterp( timeIdx ), shim_arr( :, :, pp ), PPStartStop_shim( pp, : ) ) );
    end

else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.opt_tvec = tvecOpt;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_rshp = grad_rshp;
waveform.bimag_rshp = bimag_rshp;

waveform.numVarsPerChannel_RF = numVarsPerChannel_RF;
waveform.numPPShape_RF = opt.numPPShape_RF;
waveform.orderPP_RF = opt.orderPP_RF;
waveform.numPPPeriods_RF = opt.numPPPeriods_RF;
waveform.PPIdxs_RF = opt.PPIdxs_RF;
waveform.PPStartStop_RF = opt.PPStartStop_RF;
waveform.PPVarIdxs_RF = opt.PPVarIdxs_RF;

waveform.Grad = Grad;
waveform.grad_rshp = grad_rshp;
waveform.numVarsPerChannel_grad = numVarsPerChannel_grad;
waveform.numPPShape_grad = opt.numPPShape_grad;
waveform.orderPP_grad = opt.orderPP_grad;
waveform.numPPPeriods_grad = opt.numPPPeriods_grad;
waveform.PPIdxs_grad = opt.PPIdxs_grad;
waveform.PPStartStop_grad = opt.PPStartStop_grad;
waveform.PPVarIdxs_grad = opt.PPVarIdxs_grad;

waveform.Shim = Shim;
if numZCoils > 0
    waveform.shim_rshp = shim_rshp;
    waveform.numVarsPerChannel_shim = numVarsPerChannel_shim;
    waveform.numPPShape_shim = opt.numPPShape_shim;
    waveform.orderPP_shim = opt.orderPP_shim;
    waveform.numPPPeriods_shim = opt.numPPPeriods_shim;
    waveform.PPIdxs_shim = opt.PPIdxs_shim;
    waveform.PPStartStop_shim = opt.PPStartStop_shim;
    waveform.PPVarIdxs_shim = opt.PPVarIdxs_shim;
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