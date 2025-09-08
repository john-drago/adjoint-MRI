function [ waveform ] = generatevarkTPPlotWaveform_base( p, opt )
% [ RF, Grad, Shim, RFphasor, Freq, tvec, params ] fields of waveform
% This function will generate a waveform that can be plotted based on the
% unscaled variables in pSc.

% Get unscaled parameters -- currently expecting unscaled parameters to be
% passed
% scVec = opt.scVec;
% p = pSc( : ) .* scVec;

%% Get Waveform information
numZCoils = opt.numZCoils;
numXYCoils = opt.numXYCoils;
num_kTP = opt.num_kTP;

%% Get Timing
dtperiods = p( opt.dt_idx );
tperiods = cumsum( dtperiods );

RFTiming = [ [ 0; tperiods( 2:2:end ) ], ( [ 0; tperiods( 2:2:end ) ] + dtperiods( 1:2:(end) ) ) ];
RFLength = dtperiods( 1:2:end );
RFSlewTime = opt.RFSlewTime;

blipLength = dtperiods( 2:2:end );
blipTiming = [ tperiods( 1:2:(end-2) ), ( tperiods( 1:2:(end-2) ) + dtperiods( 2:2:end ) ) ];
blipCenters = ( tperiods( 1:2:(end-2) ) + dtperiods( 2:2:end )/2 );

% Generate time vector
pulseLength = tperiods(end);
sigDigitTime = 6;

dtInterp = min( dtperiods ) / 100;
tvecInterp = round( ( 0 : dtInterp : pulseLength ), sigDigitTime, 'significant' );
dtvecInterp = dtInterp * ones( size( tvecInterp ) );
dtvecInterp( 1 ) = 0;
dtvecInterp( end ) = 0;

numTimePoints = length( tvecInterp );
dtError = eps( 1e3 );

% Determine positions of RF, grad, and slew start and stop
RF_i = zeros( num_kTP, 1, 'int32' );
RF_f = zeros( num_kTP, 1, 'int32' ); 
RF_Slew_i = zeros( num_kTP, 1, 'int32' ); 
RF_Slew_f = zeros( num_kTP, 1, 'int32' ); 

Blip_i = zeros( num_kTP-1, 1, 'int32' ); 
Blip_f = zeros( num_kTP-1, 1, 'int32' ); 

idxtol = dtError;
for nn = 1:num_kTP
    RF_i( nn ) = find( tvecInterp >= ( RFTiming(nn, 1) - idxtol ), 1, 'first' );
    RF_f( nn ) = find( tvecInterp <= ( RFTiming(nn, 2) + idxtol ), 1, 'last' );

    RF_Slew_i( nn ) = find( tvecInterp >= ( ( RFTiming( nn, 1 ) + RFSlewTime ) - idxtol  ), 1, 'first' );
    RF_Slew_f( nn ) = find( tvecInterp <= ( ( RFTiming( nn, 2 ) - RFSlewTime ) + idxtol  ), 1, 'last' );
    
    if nn ~= num_kTP
        Blip_i( nn ) = find( tvecInterp >= ( blipTiming(nn, 1) - idxtol ), 1, 'first' );
        Blip_f( nn ) = find( tvecInterp <= ( blipTiming(nn, 2) + idxtol ), 1, 'last' );
    end
end


%% Initialize Waveforms
RF = complex( zeros( numXYCoils, ( numTimePoints ) ) );
% RFphasor = complex( zeros( numXYCoils, ( numTimePoints ) ) );
Grad = zeros( 3, ( numTimePoints ) );

%% Deal with Freq
dtvec = zeros( opt.numTimePoints, 1 ); 
dtvec( opt.RF_idx ) = dtperiods( 1:2:end ) - 2*opt.RFSlewTime;
dtvec( opt.RF_Slew_i_idx ) = opt.RFSlewTime;
dtvec( opt.RF_Slew_f_idx ) = opt.RFSlewTime;
dtvec( opt.blip_idx ) = dtperiods( 2:2:end );
dtvec = transpose( dtvec(:) );

tvecOpt = cumsum( dtvec ); 
tvecOpt = ( [ 0, tvecOpt( 1:(end-1) ) ] + dtvec/2 );
tvecOpt = transpose( tvecOpt(:) );

tvecOptPad = [ 0, tvecOpt, pulseLength ];

freqInterpMethod = 'nearest';
dwxyvec_optpad = [ opt.dwxyvec(1), opt.dwxyvec, opt.dwxyvec(end) ];
Freq = interpolateWaveforms( tvecOptPad, dwxyvec_optpad / ( 2*pi ), tvecInterp, freqInterpMethod );
dwxyvec = ( 2 * pi ) * Freq;

%% Add RF points to RF array
% breal
breal_ktp_idx = opt.breal_idx;
breal_ktp = p( breal_ktp_idx );
breal_ktp_rshp = reshape( breal_ktp, [ numXYCoils, num_kTP ] );

% bimag
bimag_ktp_idx = opt.bimag_idx;
bimag_ktp = p( bimag_ktp_idx );
bimag_ktp_rshp = reshape( bimag_ktp, [ numXYCoils, num_kTP ] );

% bcomp
bcomp_ktp_rshp = complex( breal_ktp_rshp, bimag_ktp_rshp );

% Iterate over the kT-points
for kk = 1:num_kTP
    
    % rise
    rise_idx = (RF_i(kk):(RF_Slew_i( kk )-1));
    RF( :, rise_idx ) = ( interp1(...
        round([ RFTiming( kk, 1 ), (RFTiming( kk, 1 ) + RFSlewTime) ], sigDigitTime, 'significant' ),...
        [ zeros( numXYCoils, 1 ), bcomp_ktp_rshp( :, kk ) ].',...
        tvecInterp( rise_idx ) ) ).';

    % fall
    fall_idx = ((RF_Slew_f(kk)+1):(RF_f( kk )));
    RF( :, fall_idx ) = ( interp1(...
        round([ ( RFTiming( kk, 2 ) - RFSlewTime ), RFTiming( kk, 2 ) ], sigDigitTime, 'significant' ),...
        [ bcomp_ktp_rshp( :, kk ), zeros( numXYCoils, 1 ) ].',...
        tvecInterp( fall_idx ) ) ).';

    % constant
    const_idx = ((RF_Slew_i( kk )):(RF_Slew_f( kk )));
    RF( :, const_idx ) = ( repmat( bcomp_ktp_rshp( :, kk ), [ 1, length(const_idx) ] ) );
end

RFphasor = RF;

%% Add Grad points to Grad array
% grad
grad_ktp_idx = opt.grad_idx;
grad_ktp = p( grad_ktp_idx );
grad_ktp_rshp = reshape( grad_ktp, [ 3, (num_kTP-1) ] );

% Iterate over the kT-points
for kk = 1:(num_kTP-1)
    % determine interp kTP indices
    Grad_idx_interp = ( tvecInterp >= ( blipTiming(kk, 1) - dtError ) ) & ( tvecInterp <= ( blipTiming(kk, 2) + dtError ) );

    Grad( :, Grad_idx_interp ) = grad_ktp_rshp( :, kk ) * ( 1 - abs(...
        ( tvecInterp( Grad_idx_interp ) - blipCenters( kk ) ) / ( blipLength( kk ) / 2 ) ) );

end

%% Add shim points to shim array
if opt.numZCoils > 0

    Shim = zeros( numZCoils, ( numTimePoints ) );

    % shim
    shim_ktp_idx = opt.shim_idx;
    shim_ktp = p( shim_ktp_idx );
    shim_ktp_rshp = reshape( shim_ktp, [ numZCoils, (num_kTP-1) ] );

    % Iterate over the kT-points
    for kk = 1:(num_kTP-1)

        % determine interp kTP indices
        Shim_idx_interp = ( tvecInterp >= ( blipTiming(kk, 1) - dtError ) ) &...
            ( tvecInterp <= ( blipTiming(kk, 2) + dtError ) );

        Shim( :, Shim_idx_interp ) = shim_ktp_rshp( :, kk ) * ( 1 - abs(...
            ( tvecInterp( Shim_idx_interp ) - blipCenters( kk ) ) / ( blipLength( kk ) / 2 ) ) );

    end

else

    Shim = zeros( 1, ( numTimePoints ) );

end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_ktp_rshp = breal_ktp_rshp;
waveform.bimag_ktp_rshp = bimag_ktp_rshp;
waveform.bcomp_ktp_rshp = bcomp_ktp_rshp;

waveform.Grad = Grad;
waveform.grad_ktp_rshp = grad_ktp_rshp;

waveform.Shim = Shim;
if opt.numZCoils > 0
    waveform.shim_ktp_rshp = shim_ktp_rshp;
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

waveform.num_kTP = num_kTP;
waveform.blipLength = blipLength;
waveform.blipTiming = blipTiming;
waveform.blipCenters = blipCenters;
waveform.RFLength = RFLength;

waveform.Z0 = opt.Z0;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %