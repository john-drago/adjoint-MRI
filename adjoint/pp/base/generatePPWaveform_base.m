function [ waveform ] = generatePPWaveform_base( p, opt )
% [ RF, Grad, Shim, RFphasor, Freq, tvec, params ] fields of waveform
% Generate waveform for the hard pulse integration.

% Get unscaled parameters -- currently expecting unscaled parameters to be
% passed
% scVec = opt.scVec;
% p = pSc( : ) .* scVec;
p = p( : );

% Initialize Waveforms
numZCoils = opt.numZCoils;
numXYCoils = opt.numXYCoils;
numTimePoints = opt.numTimePoints;
tvec = opt.tvec;
dtvec = opt.dtvec;
dwxyvec = opt.dwxyvec;

% brealphasor = zeros( numXYCoils, numTimePoints );
% bimagphasor = zeros( numXYCoils, numTimePoints );
% breal = zeros( numXYCoils, numTimePoints );
% bimag = zeros( numXYCoils, numTimePoints );
% Grad = zeros( 3, numTimePoints );

%% Add phasor points for RF
numVarsPerChannel_RF = opt.numVarsPerChannel_RF;

% breal
breal_idx = opt.breal_idx;
breal_vec = p( breal_idx );
breal_rshp = reshape( breal_vec, [ numVarsPerChannel_RF, numXYCoils ] );

% bimag
bimag_idx = opt.bimag_idx;
bimag_vec = p( bimag_idx );
bimag_rshp = reshape( bimag_vec, [ numVarsPerChannel_RF, numXYCoils ] );

breal = transpose( opt.varsToTimepoints_RF * breal_rshp );
bimag = transpose( opt.varsToTimepoints_RF * bimag_rshp );

brealphasor = breal;
bimagphasor = bimag;

%% Add Grad points
numVarsPerChannel_grad = opt.numVarsPerChannel_grad;

% grad
grad_idx = opt.grad_idx;
grad_vec = p( grad_idx );
grad_rshp = reshape( grad_vec, [ numVarsPerChannel_grad, 3 ] );

Grad = transpose( opt.varsToTimepoints_grad * grad_rshp );

%% Add Shim points
if numZCoils > 0

    numVarsPerChannel_shim = opt.numVarsPerChannel_shim;

    % shim
    shim_idx = opt.shim_idx;
    shim_vec = p( shim_idx );
    shim_rshp = reshape( shim_vec, [ numVarsPerChannel_shim, numZCoils ] );

    Shim = transpose( opt.varsToTimepoints_shim * shim_rshp );

else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.breal_rshp = breal_rshp;
waveform.bimag_rshp = bimag_rshp;

waveform.numVarsPerChannel_RF = numVarsPerChannel_RF;
waveform.numPPShape_RF = opt.numPPShape_RF;
waveform.orderPP_RF = opt.orderPP_RF;
waveform.numPPPeriods_RF = opt.numPPPeriods_RF;
waveform.PPIdxs_RF = opt.PPIdxs_RF;
waveform.PPStartStop_RF = opt.PPStartStop_RF;
waveform.PPVarIdxs_RF = opt.PPVarIdxs_RF;
waveform.varsToTimepoints_RF = opt.varsToTimepoints_RF;

waveform.grad_rshp = grad_rshp;
waveform.numVarsPerChannel_grad = numVarsPerChannel_grad;
waveform.numPPShape_grad = opt.numPPShape_grad;
waveform.orderPP_grad = opt.orderPP_grad;
waveform.numPPPeriods_grad = opt.numPPPeriods_grad;
waveform.PPIdxs_grad = opt.PPIdxs_grad;
waveform.PPStartStop_grad = opt.PPStartStop_grad;
waveform.PPVarIdxs_grad = opt.PPVarIdxs_grad;
waveform.varsToTimepoints_grad = opt.varsToTimepoints_grad;

if numZCoils > 0
    waveform.shim_rshp = shim_rshp;
    waveform.numVarsPerChannel_shim = numVarsPerChannel_shim;
    waveform.numPPShape_shim = opt.numPPShape_shim;
    waveform.orderPP_shim = opt.orderPP_shim;
    waveform.numPPPeriods_shim = opt.numPPPeriods_shim;
    waveform.PPIdxs_shim = opt.PPIdxs_shim;
    waveform.PPStartStop_shim = opt.PPStartStop_shim;
    waveform.PPVarIdxs_shim = opt.PPVarIdxs_shim;
    waveform.varsToTimepoints_shim = opt.varsToTimepoints_shim;
end

% if matches( opt.structtype, 'opt' )
%     if opt.useGPU
%         waveform.shapeFnValsTimePoints_RF_rep = opt.gpu_shapeFnValsTimePoints_RF_rep;
%         waveform.shapeFnValsTimePoints_grad_rep = opt.gpu_shapeFnValsTimePoints_grad_rep;
%         if numZCoils > 0
%             waveform.shapeFnValsTimePoints_shim_rep = opt.gpu_shapeFnValsTimePoints_shim_rep;
%         end
%     else
%         waveform.shapeFnValsTimePoints_RF_rep = opt.shapeFnValsTimePoints_RF_rep;
%         waveform.shapeFnValsTimePoints_grad_rep = opt.shapeFnValsTimePoints_grad_rep;
%         if numZCoils > 0
%             waveform.shapeFnValsTimePoints_shim_rep = opt.shapeFnValsTimePoints_shim_rep;
%         end
%     end
% end

end