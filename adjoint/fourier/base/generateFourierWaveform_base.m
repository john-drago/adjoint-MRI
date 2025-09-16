function [ waveform ] = generateFourierWaveform_base( p, opt )
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
omega0 = opt.omega0;

% brealphasor = zeros( numXYCoils, numTimePoints );
% bimagphasor = zeros( numXYCoils, numTimePoints );
% breal = zeros( numXYCoils, numTimePoints );
% bimag = zeros( numXYCoils, numTimePoints );
% Grad = zeros( 3, numTimePoints );

%% Add phasor points for RF
numFourier_RF = opt.numFourier_RF;

% breal
breal_idx = opt.breal_idx;
breal_vec = p( breal_idx );
breal_rshp = reshape( breal_vec, [ numFourier_RF, numXYCoils ] );

% bimag
bimag_idx = opt.bimag_idx;
bimag_vec = p( bimag_idx );
bimag_rshp = reshape( bimag_vec, [ numFourier_RF, numXYCoils ] );

bcomp_rshp = complex(...
    breal_rshp, bimag_rshp );

if strcmpi( opt.structtype, "opt" )
    
    bcomp = opt.FBRF * bcomp_rshp;
    breal = real( bcomp );
    bimag = imag( bcomp );

elseif strcmpi( opt.structtype, "val" )

    FBRF = getFourierRFBasisMatrix( opt.orderFourier_RF, tvec, omega0 );

    bcomp = FBRF * bcomp_rshp;
    breal = real( bcomp );
    bimag = imag( bcomp );
    
end

breal = transpose( breal );
bimag = transpose( bimag );

brealphasor = breal;
bimagphasor = bimag;

%% Add Grad points
numFourier_grad = opt.numFourier_grad;

% grad
grad_idx = opt.grad_idx;
grad_vec = p( grad_idx );
grad_rshp = reshape( grad_vec, [ numFourier_grad, 3 ] );

if strcmpi( opt.structtype, "opt" )
    Grad = opt.FBgrad * grad_rshp;
elseif strcmpi( opt.structtype, "val" )

    FBgrad = getFourierBasisMatrix( opt.orderFourier_grad, tvec, omega0 );

    Grad = FBgrad * grad_rshp;

end

Grad = transpose( Grad );

%% Add Shim points
if numZCoils > 0
    
    numFourier_shim = opt.numFourier_shim;
    shim_idx = opt.shim_idx;
    shim_vec = p( shim_idx );
    shim_rshp = reshape( shim_vec, [ numFourier_shim, numZCoils ] );

    if strcmpi( opt.structtype, "opt" )
        Shim = opt.FBshim * shim_rshp;
    elseif strcmpi( opt.structtype, "val" )

        FBshim = getFourierBasisMatrix( opt.orderFourier_shim, tvec, omega0 );

        Shim = FBshim * shim_rshp;
    end

    Shim = transpose( Shim );
    
else
    Shim = zeros( 1, ( numTimePoints ) );
end



%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.bcomp_rshp = bcomp_rshp;
waveform.numFourier_RF = numFourier_RF;

waveform.grad_rshp = grad_rshp;
waveform.numFourier_grad = numFourier_grad;

if numZCoils > 0
    waveform.shim_rshp = shim_rshp;
    waveform.numFourier_shim = numFourier_shim;
end

if matches( opt.structtype, 'opt' )
    if opt.useGPU
        waveform.FBRF = opt.gpu_FBRF;
        waveform.FBgrad = opt.gpu_FBgrad;
        if opt.numZCoils > 0
            waveform.FBshim = opt.gpu_FBshim;
        end
    else
        waveform.FBRF = opt.FBRF;
        waveform.FBgrad = opt.FBgrad;
        if opt.numZCoils > 0
            waveform.FBshim = opt.FBshim;
        end
    end
end

end