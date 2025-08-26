function [ waveform ] = generatekTPWaveform_base( p, opt )
% [ RF, Grad, Shim, RFphasor, Freq, tvec, params ] fields of waveform
% Generate waveform for the hard pulse integration.

% Get unscaled parameters -- currently expecting unscaled parameters to be
% passed
% scVec = opt.scVec;
% p = pSc( : ) .* scVec;
p = p( : );

%% Get Timing
dtvec = opt.dtvec;
tvec = opt.tvec;
dwxyvec = opt.dwxyvec;

RFTiming = opt.RFTiming;
RFLength = opt.RFLength;
RFSlewTime = opt.RFSlewTime;
blipLength = opt.blipLength;
blipTiming = opt.blipTiming;
blipCenters = opt.blipCenters;

%% Initialize Waveforms
numZCoils = opt.numZCoils;
numXYCoils = opt.numXYCoils;
num_kTP = opt.num_kTP;
numTimePoints = opt.numTimePoints;

% brealphasor = zeros( numXYCoils, numTimePoints );
% bimagphasor = zeros( numXYCoils, numTimePoints );
breal = zeros( numXYCoils, numTimePoints );
bimag = zeros( numXYCoils, numTimePoints );
Grad = zeros( 3, numTimePoints );

%% Add phasor points for RF
% breal
breal_ktp_idx = opt.breal_idx;
breal_ktp = p( breal_ktp_idx );
breal_ktp_rshp = reshape( breal_ktp, [ numXYCoils, num_kTP ] );

% bimag
bimag_ktp_idx = opt.bimag_idx;
bimag_ktp = p( bimag_ktp_idx );
bimag_ktp_rshp = reshape( bimag_ktp, [ numXYCoils, num_kTP ] );

% Iterate over the kT-points

if strcmpi( opt.structtype, "opt" )

    breal( :, opt.RF_idx ) = breal_ktp_rshp;
    bimag( :, opt.RF_idx ) = bimag_ktp_rshp;

    breal( :, opt.RF_Slew_i_idx ) = 0.5 * breal_ktp_rshp;
    bimag( :, opt.RF_Slew_i_idx ) = 0.5 * bimag_ktp_rshp;

    breal( :, opt.RF_Slew_f_idx ) = 0.5 * breal_ktp_rshp;
    bimag( :, opt.RF_Slew_f_idx ) = 0.5 * bimag_ktp_rshp;

elseif strcmpi( opt.structtype, "val" )

    breal( :, opt.RF_idx ) = breal_ktp_rshp;
    bimag( :, opt.RF_idx ) = bimag_ktp_rshp;

    breal( :, opt.RF_Slew_i_idx ) = 0.5 * breal_ktp_rshp;
    bimag( :, opt.RF_Slew_i_idx ) = 0.5 * bimag_ktp_rshp;

    breal( :, opt.RF_Slew_f_idx ) = 0.5 * breal_ktp_rshp;
    bimag( :, opt.RF_Slew_f_idx ) = 0.5 * bimag_ktp_rshp;

end

brealphasor = breal;
bimagphasor = bimag;

%% Add Grad points
% grad
grad_ktp_idx = opt.grad_idx;
grad_ktp = p( grad_ktp_idx );
grad_ktp_rshp = reshape( grad_ktp, [ 3, (num_kTP-1) ] );

% Iterate over the kT-points
if strcmpi( opt.structtype, "opt" )

    Grad( :, opt.blip_idx ) = 0.5 * grad_ktp_rshp;

elseif strcmpi( opt.structtype, "val" )

    Grad( :, opt.blip_idx ) = 0.5 * grad_ktp_rshp;
end

%% Add Shim points
if opt.numZCoils > 0
    Shim = zeros( numZCoils, numTimePoints );

    % shim
    shim_ktp_idx = opt.shim_idx;
    shim_ktp = p( shim_ktp_idx );
    shim_ktp_rshp = reshape( shim_ktp, [ numZCoils, (num_kTP-1) ] );

    % Iterate over the kT-points
    if strcmpi( opt.structtype, "opt" )

        Shim( :, opt.blip_idx ) = 0.5 * shim_ktp_rshp;

    elseif strcmpi( opt.structtype, "val" )

        Shim( :, opt.blip_idx ) = 0.5 * shim_ktp_rshp;
    end
else
    Shim = zeros( 1, numTimePoints );
end

%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.breal_ktp_rshp = breal_ktp_rshp;
waveform.bimag_ktp_rshp = bimag_ktp_rshp;

waveform.grad_ktp_rshp = grad_ktp_rshp;

if opt.numZCoils > 0
    waveform.shim_ktp_rshp = shim_ktp_rshp;
    waveform.shimSlewTime = opt.shimSlewTime;
end

waveform.num_kTP = num_kTP;
waveform.RFTiming = RFTiming;
waveform.RFLength = RFLength;
waveform.RFSlewTime = RFSlewTime;
waveform.blipLength = blipLength;
waveform.blipTiming = blipTiming;
waveform.blipCenters = blipCenters;

waveform.gradSlewTime = opt.gradSlewTime;

end