function [ waveform ] = generatePWCWaveform_base( p, opt )
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
breal = zeros( numXYCoils, numTimePoints );
bimag = zeros( numXYCoils, numTimePoints );
Grad = zeros( 3, numTimePoints );

%% Add phasor points for RF
% breal
breal_idx = opt.breal_idx;
breal_vec = p( breal_idx );
breal_rshp = reshape( breal_vec, [ numTimePoints, numXYCoils ] ).';

% bimag
bimag_idx = opt.bimag_idx;
bimag_vec = p( bimag_idx );
bimag_rshp = reshape( bimag_vec, [ numTimePoints, numXYCoils ] ).';

if strcmpi( opt.structtype, "opt" )
    
    breal( :, 1:end ) = breal_rshp;
    bimag( :, 1:end ) = bimag_rshp;

elseif strcmpi( opt.structtype, "val" )
    
    breal( :, 1:end ) = interp1( [ 0, opt.opt_tvec, opt.pulseLength ].', [ zeros( numXYCoils, 1), breal_rshp, zeros( numXYCoils, 1) ].', opt.tvec.' ).';
    bimag( :, 1:end ) = interp1( [ 0, opt.opt_tvec, opt.pulseLength ].', [ zeros( numXYCoils, 1), bimag_rshp, zeros( numXYCoils, 1) ].', opt.tvec.' ).';
    
end

brealphasor = breal;
bimagphasor = bimag;

%% Add Grad points
% grad
grad_idx = opt.grad_idx;
grad_vec = p( grad_idx );
grad_rshp = transpose( reshape( grad_vec, [ numTimePoints, 3 ] ) );

if strcmpi( opt.structtype, "opt" )
    Grad( :, 1:end ) = grad_rshp;
elseif strcmpi( opt.structtype, "val" )
    Grad( :, 1:end ) = transpose( interp1( transpose( [ 0, opt.opt_tvec, opt.pulseLength ] ), transpose( [ zeros( 3, 1), grad_rshp, zeros( 3, 1) ] ), transpose( opt.tvec ) ) );
end

%% Add Shim points
if numZCoils > 0
    Shim = zeros( numZCoils, numTimePoints );

    shim_idx = opt.shim_idx;
    shim_vec = p( shim_idx );

    shim_rshp = reshape( shim_vec, [ numTimePoints, numZCoils ] ).';

    if strcmpi( opt.structtype, "opt" )
        Shim( :, 1:end ) = shim_rshp;
    elseif strcmpi( opt.structtype, "val" )
        Shim( :, 1:end ) = transpose( interp1( transpose( [ 0, opt.opt_tvec, opt.pulseLength ] ), transpose( [ zeros( numZCoils, 1), shim_rshp, zeros( numZCoils, 1) ] ), transpose( opt.tvec ) ) );
    end
else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.breal_vec = breal_vec;
waveform.bimag_vec = bimag_vec;

waveform.grad_vec = grad_vec;

if numZCoils > 0
    waveform.shim_vec = shim_vec;
end

end