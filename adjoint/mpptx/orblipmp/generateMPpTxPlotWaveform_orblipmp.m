function [ waveform ] = generateMPpTxPlotWaveform_orblipmp( p, opt )
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
dtInterp = opt.dt/2;
pulseLength = opt.pulseLength;
sigDigitTime = 6;
tvecInterp = round( 0 : dtInterp : pulseLength, sigDigitTime, 'significant' );
dtvecInterp = dtInterp * ones( size( tvecInterp ) );
dtvecInterp( 1 ) = 0;
dtvecInterp( end ) = 0;

tvecOrig = opt.tvec;

numTimePoints = length( tvecInterp );
dtError = eps( 1e3 );

RF = complex( zeros( numXYCoils, ( numTimePoints ) ) );
RFphasor = complex( zeros( numXYCoils, ( numTimePoints ) ) );
Freq = zeros( 1, ( numTimePoints ) );
Grad = zeros( 3, ( numTimePoints ) );

% Get constants
dwxy_mpptx = opt.dwxy_mpptx;
wz_mpptx = opt.wz_mpptx;
dfxy_mpptx = opt.dfxy_mpptx;
fz_mpptx = opt.fz_mpptx;

%% Add RF points to RF array
breal_ORSP_p_idx = opt.breal_ORSP_idx;
bimag_ORSP_p_idx = opt.bimag_ORSP_idx;
breal_MPSP_p_idx = opt.breal_MPSP_idx;
bimag_MPSP_p_idx = opt.bimag_MPSP_idx;

breal_ORSP_vec = p( breal_ORSP_p_idx );
bimag_ORSP_vec = p( bimag_ORSP_p_idx );
breal_MPSP_vec = p( breal_MPSP_p_idx );
bimag_MPSP_vec = p( bimag_MPSP_p_idx );

% Add points during the on-resonance subpulse
ORSP_idxs = opt.ORSP_i : opt.ORSP_f;

breal_ORSP_rshp = reshape( breal_ORSP_vec, [ length( ORSP_idxs ), numXYCoils ] ).';
bimag_ORSP_rshp = reshape( bimag_ORSP_vec, [ length( ORSP_idxs ), numXYCoils ] ).';

ORSP_idx_interp = ( tvecInterp >= ( opt.tStORSP - dtError ) & ( tvecInterp <= ( opt.tEndORSP + dtError )  ) );
tvecInterp_ORSP = tvecInterp( ORSP_idx_interp );
breal_ORSP_interp = interp1( [0, tvecOrig( ORSP_idxs ), opt.tEndORSP], [ zeros( numXYCoils, 1), breal_ORSP_rshp, zeros( numXYCoils, 1) ].', tvecInterp_ORSP ).';
bimag_ORSP_interp = interp1( [0, tvecOrig( ORSP_idxs ), opt.tEndORSP], [ zeros( numXYCoils, 1), bimag_ORSP_rshp, zeros( numXYCoils, 1) ].', tvecInterp_ORSP ).';

RF( :, ORSP_idx_interp ) = breal_ORSP_interp + 1j * bimag_ORSP_interp;
RFphasor( :, ORSP_idx_interp ) = breal_ORSP_interp + 1j * bimag_ORSP_interp;

% Add points during the multiphoton subpulse
MPSP_idxs = opt.MPSP_i : opt.MPSP_f;

breal_MPSP_rshp = reshape( breal_MPSP_vec, [ length( MPSP_idxs ), numXYCoils ] ).';
bimag_MPSP_rshp = reshape( bimag_MPSP_vec, [ length( MPSP_idxs ), numXYCoils ] ).';

MPSP_idx_interp = ( tvecInterp >= ( opt.tStMPSP - dtError ) & ( tvecInterp <= ( opt.tEndMPSP + dtError )  ) );
tvecInterp_MPSP = tvecInterp( MPSP_idx_interp );
breal_MPSP_interp = interp1( [opt.tStMPSP, tvecOrig( MPSP_idxs ), opt.tEndMPSP], [ zeros( numXYCoils, 1), breal_MPSP_rshp, zeros( numXYCoils, 1) ].', tvecInterp_MPSP ).';
bimag_MPSP_interp = interp1( [opt.tStMPSP, tvecOrig( MPSP_idxs ), opt.tEndMPSP], [ zeros( numXYCoils, 1), bimag_MPSP_rshp, zeros( numXYCoils, 1) ].', tvecInterp_MPSP ).';

RFphasor( :, MPSP_idx_interp ) = breal_MPSP_interp + 1j * bimag_MPSP_interp;
RF( :, MPSP_idx_interp ) = RFphasor( :, MPSP_idx_interp ) .* exp( 1j * dwxy_mpptx * ( tvecInterp_MPSP - opt.tStMPSP ) );

Freq( :, MPSP_idx_interp ) = dfxy_mpptx;

%% Add Grad points to Grad array

% Add points during the multiphoton subpulse
Gx_ORSP_p_idx = opt.Gx_ORSP_idx;
Gy_ORSP_p_idx = opt.Gy_ORSP_idx;
Gz_ORSP_p_idx = opt.Gz_ORSP_idx;

Gx_ORSP = p( Gx_ORSP_p_idx );
Gy_ORSP = p( Gy_ORSP_p_idx );
Gz_ORSP = p( Gz_ORSP_p_idx );

Grad( 1, ORSP_idx_interp ) = interp1( [0, tvecOrig( ORSP_idxs ), opt.tEndORSP], [ 0, Gx_ORSP.', 0 ], tvecInterp_ORSP );
Grad( 2, ORSP_idx_interp ) = interp1( [0, tvecOrig( ORSP_idxs ), opt.tEndORSP], [ 0, Gy_ORSP.', 0 ], tvecInterp_ORSP );
Grad( 3, ORSP_idx_interp ) = interp1( [0, tvecOrig( ORSP_idxs ), opt.tEndORSP], [ 0, Gz_ORSP.', 0 ], tvecInterp_ORSP );

% Add points during the blip period
Gx_Blip_p_idx = opt.Gx_Blip_idx;
Gy_Blip_p_idx = opt.Gy_Blip_idx;
Gz_Blip_p_idx = opt.Gz_Blip_idx;

Gx_Blip = p( Gx_Blip_p_idx );
Gy_Blip = p( Gy_Blip_p_idx );
Gz_Blip = p( Gz_Blip_p_idx );

G_Blip = [ Gx_Blip; Gy_Blip; Gz_Blip ];

G_Blip_slew_i_idxs = opt.Blip_i:( opt.Blip_Grad_Slew_i-1 );
G_Blip_slew_f_idxs = ( opt.Blip_Grad_Slew_f+1 ):( opt.Blip_f );
G_Blip_const_idxs = ( opt.Blip_Grad_Slew_i ) : opt.Blip_Grad_Slew_f;

Blip_Grad_slew_i_idx_interp = ( tvecInterp >= ( opt.tStBlip - dtError ) & ( tvecInterp <= ( ( opt.tStBlip + opt.BlipGradSlewTime ) + dtError )  ) );
tvecInterp_Grad_Blip_slew_i = tvecInterp( Blip_Grad_slew_i_idx_interp );
Blip_Grad_slew_f_idx_interp = ( tvecInterp >= ( ( opt.tEndBlip - opt.BlipGradSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndBlip ) + dtError )  ) );
tvecInterp_Grad_Blip_slew_f = tvecInterp( Blip_Grad_slew_f_idx_interp );
Blip_Grad_const_idx_interp = ( tvecInterp >= ( ( opt.tStBlip + opt.BlipGradSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndBlip - opt.BlipGradSlewTime ) + dtError )  ) );
tvecInterp_Blip_Grad_const = tvecInterp( Blip_Grad_const_idx_interp );

% slew waveforms
Grad( :, Blip_Grad_slew_i_idx_interp ) = transpose( interp1(...
    round( [ opt.tStBlip, tvecOrig( G_Blip_slew_i_idxs ), (opt.tStBlip + opt.BlipGradSlewTime) ], sigDigitTime, 'significant'),...
    [ zeros( 3, 1 ), ( G_Blip / opt.BlipGradSlewTime ) * ( tvecOrig( G_Blip_slew_i_idxs ) - opt.tStBlip ), G_Blip ].',...
    tvecInterp_Grad_Blip_slew_i) );

Grad( :, Blip_Grad_slew_f_idx_interp ) = transpose( interp1(...
    round( [ ( opt.tEndBlip - opt.BlipGradSlewTime ), tvecOrig( G_Blip_slew_f_idxs ), (opt.tEndBlip) ], sigDigitTime, 'significant' ),...
    [ G_Blip, -( G_Blip / opt.BlipGradSlewTime ) * ( tvecOrig( G_Blip_slew_f_idxs ) - opt.tEndBlip ), zeros( 3, 1 ) ].',...
    tvecInterp_Grad_Blip_slew_f) );

% constant blip waveforms
Grad( :, Blip_Grad_const_idx_interp ) = transpose( interp1(...
    round( [ ( opt.tStBlip + opt.BlipGradSlewTime ), tvecOrig( G_Blip_const_idxs ), (opt.tEndBlip - opt.BlipGradSlewTime) ], sigDigitTime, 'significant' ),...
    repmat( G_Blip, [ 1, length( tvecOrig( G_Blip_const_idxs ) ) + 2 ] ).',...
    tvecInterp_Blip_Grad_const) );

% Add points during the MPSP
Gxreal_MPSP_p_idx = opt.Gxreal_MPSP_idx;
Gximag_MPSP_p_idx = opt.Gximag_MPSP_idx;

Gyreal_MPSP_p_idx = opt.Gyreal_MPSP_idx;
Gyimag_MPSP_p_idx = opt.Gyimag_MPSP_idx;

Gzreal_MPSP_p_idx = opt.Gzreal_MPSP_idx;
Gzimag_MPSP_p_idx = opt.Gzimag_MPSP_idx;

Gx_real_MPSP = p( Gxreal_MPSP_p_idx );
Gx_imag_MPSP = p( Gximag_MPSP_p_idx );
Gy_real_MPSP = p( Gyreal_MPSP_p_idx );
Gy_imag_MPSP = p( Gyimag_MPSP_p_idx );
Gz_real_MPSP = p( Gzreal_MPSP_p_idx );
Gz_imag_MPSP = p( Gzimag_MPSP_p_idx );

G_real_MPSP = [ Gx_real_MPSP; Gy_real_MPSP; Gz_real_MPSP ];
G_imag_MPSP = [ Gx_imag_MPSP; Gy_imag_MPSP; Gz_imag_MPSP ];
G_mag_MPSP = sqrt( G_real_MPSP.^2 + G_imag_MPSP.^2 );
G_ph_MPSP = atan2( G_imag_MPSP, G_real_MPSP );

G_MPSP_slew_i_idxs = opt.MPSP_i:( opt.MPSP_Grad_Slew_i-1 );
G_MPSP_slew_f_idxs = ( opt.MPSP_Grad_Slew_f+1 ):( opt.MPSP_f );
G_MPSP_wave_idxs = ( opt.MPSP_Grad_Slew_i ) : opt.MPSP_Grad_Slew_f;

MPSP_Grad_slew_i_idx_interp = ( tvecInterp >= ( opt.tStMPSP - dtError ) & ( tvecInterp <= ( ( opt.tStMPSP + opt.MPSPGradSlewTime ) + dtError )  ) );
tvecInterp_Grad_MPSP_slew_i = tvecInterp( MPSP_Grad_slew_i_idx_interp );
MPSP_Grad_slew_f_idx_interp = ( tvecInterp >= ( ( opt.tEndMPSP - opt.MPSPGradSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndMPSP ) + dtError )  ) );
tvecInterp_Grad_MPSP_slew_f = tvecInterp( MPSP_Grad_slew_f_idx_interp );
MPSP_Grad_const_idx_interp = ( tvecInterp >= ( ( opt.tStMPSP + opt.MPSPGradSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndMPSP - opt.MPSPGradSlewTime ) + dtError )  ) );
tvecInterp_MPSP_Grad_const = tvecInterp( MPSP_Grad_const_idx_interp );

% slew waveforms
Grad( :, MPSP_Grad_slew_i_idx_interp ) = ( interp1(...
    round( [ opt.tStMPSP, tvecOrig( G_MPSP_slew_i_idxs ), (opt.tStMPSP + opt.MPSPGradSlewTime) ], sigDigitTime, 'significant' ),...
    [ zeros( 3, 1 ), ( ( G_mag_MPSP .* cos( wz_mpptx * ( opt.MPSPGradSlewTime ) + G_ph_MPSP ) ) / opt.MPSPGradSlewTime ) * ( tvecOrig( G_MPSP_slew_i_idxs ) - opt.tStMPSP ), ( G_mag_MPSP .* cos( wz_mpptx * ( opt.MPSPGradSlewTime ) + G_ph_MPSP ) ) ].',...
    tvecInterp_Grad_MPSP_slew_i) ).';

Grad( :, MPSP_Grad_slew_f_idx_interp ) = ( interp1(...
    round( [ ( opt.tEndMPSP - opt.MPSPGradSlewTime ), tvecOrig( G_MPSP_slew_f_idxs ), (opt.tEndMPSP) ], sigDigitTime, 'significant' ),...
    [ ( G_mag_MPSP .* cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPGradSlewTime ) + G_ph_MPSP ) ), -( ( G_mag_MPSP .* cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPGradSlewTime ) + G_ph_MPSP ) ) / opt.MPSPGradSlewTime ) * ( tvecOrig( G_MPSP_slew_f_idxs ) - opt.tEndMPSP ), zeros( 3, 1 ) ].',...
    tvecInterp_Grad_MPSP_slew_f) ).';

% constant MPSP waveforms
Grad( :, MPSP_Grad_const_idx_interp ) = ( interp1(...
    round( [ ( opt.tStMPSP + opt.MPSPGradSlewTime ), tvecOrig( G_MPSP_wave_idxs ), (opt.tEndMPSP - opt.MPSPGradSlewTime) ], sigDigitTime, 'significant' ),...
    (G_mag_MPSP .* cos( wz_mpptx * (  [ ( opt.tStMPSP + opt.MPSPGradSlewTime ), tvecOrig( G_MPSP_wave_idxs ), (opt.tEndMPSP - opt.MPSPGradSlewTime) ] - opt.tStMPSP ) + G_ph_MPSP )).',...
    tvecInterp_MPSP_Grad_const) ).';

%% Add shim points to shim array
if numZCoils > 0

    Shim = zeros( numZCoils, ( numTimePoints ) );

    % Add points during the blip period
    Shim_Blip_p_idx = opt.shim_Blip_idx;

    Shim_Blip = p( Shim_Blip_p_idx );

    Shim_Blip_slew_i_idxs = opt.Blip_i:( opt.Blip_Shim_Slew_i-1 );
    Shim_Blip_slew_f_idxs = ( opt.Blip_Shim_Slew_f+1 ):( opt.Blip_f );
    Shim_Blip_const_idxs = ( opt.Blip_Shim_Slew_i ) : opt.Blip_Shim_Slew_f;

    Blip_Shim_slew_i_idx_interp = ( tvecInterp >= ( opt.tStBlip - dtError ) & ( tvecInterp <= ( ( opt.tStBlip + opt.BlipShimSlewTime ) + dtError )  ) );
    tvecInterp_Shim_Blip_slew_i = tvecInterp( Blip_Shim_slew_i_idx_interp );
    Blip_Shim_slew_f_idx_interp = ( tvecInterp >= ( ( opt.tEndBlip - opt.BlipShimSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndBlip ) + dtError )  ) );
    tvecInterp_Shim_Blip_slew_f = tvecInterp( Blip_Shim_slew_f_idx_interp );
    Blip_Shim_const_idx_interp = ( tvecInterp >= ( ( opt.tStBlip + opt.BlipShimSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndBlip - opt.BlipShimSlewTime ) + dtError )  ) );
    tvecInterp_Blip_Shim_const = tvecInterp( Blip_Shim_const_idx_interp );

    % slew waveforms
    Shim( :, Blip_Shim_slew_i_idx_interp ) = ( interp1(...
        round([ opt.tStBlip, tvecOrig( Shim_Blip_slew_i_idxs ), (opt.tStBlip + opt.BlipShimSlewTime) ], sigDigitTime, 'significant' ),...
        [ zeros( opt.numZCoils, 1 ), ( Shim_Blip / opt.BlipShimSlewTime ) * ( tvecOrig( Shim_Blip_slew_i_idxs ) - opt.tStBlip ), Shim_Blip ].',...
        tvecInterp_Shim_Blip_slew_i) ).';

    Shim( :, Blip_Shim_slew_f_idx_interp ) = ( interp1(...
        round([ ( opt.tEndBlip - opt.BlipShimSlewTime ), tvecOrig( Shim_Blip_slew_f_idxs ), (opt.tEndBlip) ], sigDigitTime, 'significant' ),...
        [ Shim_Blip, -( Shim_Blip / opt.BlipShimSlewTime ) * ( tvecOrig( Shim_Blip_slew_f_idxs ) - opt.tEndBlip ), zeros( opt.numZCoils, 1 ) ].',...
        tvecInterp_Shim_Blip_slew_f) ).';

    % constant blip waveforms
    Shim( :, Blip_Shim_const_idx_interp ) = ( interp1(...
        round([ ( opt.tStBlip + opt.BlipShimSlewTime ), tvecOrig( Shim_Blip_const_idxs ), (opt.tEndBlip - opt.BlipShimSlewTime) ], sigDigitTime, 'significant' ),...
        repmat( Shim_Blip, [ 1, length( tvecOrig( Shim_Blip_const_idxs ) ) + 2 ] ).',...
        tvecInterp_Blip_Shim_const) ).';

    % Add points during the MPSP
    Shimreal_MPSP_p_idx = opt.shimreal_MPSP_idx;
    Shimimag_MPSP_p_idx = opt.shimimag_MPSP_idx;

    Shimreal_MPSP = p( Shimreal_MPSP_p_idx );
    Shimimag_MPSP = p( Shimimag_MPSP_p_idx );

    Shim_mag_MPSP = sqrt( Shimreal_MPSP.^2 + Shimimag_MPSP.^2 );
    Shim_ph_MPSP = atan2( Shimimag_MPSP, Shimreal_MPSP );

    Shim_MPSP_slew_i_idxs = opt.MPSP_i:( opt.MPSP_Shim_Slew_i-1 );
    Shim_MPSP_slew_f_idxs = ( opt.MPSP_Shim_Slew_f+1 ):( opt.MPSP_f );
    Shim_MPSP_wave_idxs = ( opt.MPSP_Shim_Slew_i ) : opt.MPSP_Shim_Slew_f;

    MPSP_Shim_slew_i_idx_interp = ( tvecInterp >= ( opt.tStMPSP ) & ( tvecInterp <= ( ( opt.tStMPSP + opt.MPSPShimSlewTime ) + dtError )  ) );
    tvecInterp_Shim_MPSP_slew_i = tvecInterp( MPSP_Shim_slew_i_idx_interp );
    MPSP_Shim_slew_f_idx_interp = ( tvecInterp >= ( ( opt.tEndMPSP - opt.MPSPShimSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndMPSP ) + dtError )  ) );
    tvecInterp_Shim_MPSP_slew_f = tvecInterp( MPSP_Shim_slew_f_idx_interp );
    MPSP_Shim_const_idx_interp = ( tvecInterp >= ( ( opt.tStMPSP + opt.MPSPShimSlewTime ) - dtError ) & ( tvecInterp <= ( ( opt.tEndMPSP - opt.MPSPShimSlewTime ) + dtError )  ) );
    tvecInterp_MPSP_Shim_const = tvecInterp( MPSP_Shim_const_idx_interp );

    % slew waveforms
    Shim( :, MPSP_Shim_slew_i_idx_interp ) = ( interp1(...
        round( [ opt.tStMPSP, tvecOrig( Shim_MPSP_slew_i_idxs ), (opt.tStMPSP + opt.MPSPShimSlewTime) ], sigDigitTime, 'significant' ),...
        [ zeros( opt.numZCoils, 1 ),...
        ( ( Shim_mag_MPSP .* cos( wz_mpptx * ( opt.MPSPShimSlewTime ) + Shim_ph_MPSP ) ) / opt.MPSPShimSlewTime ) * ( tvecOrig( Shim_MPSP_slew_i_idxs ) - opt.tStMPSP ),...
        ( Shim_mag_MPSP .* cos( wz_mpptx * ( opt.MPSPShimSlewTime ) + Shim_ph_MPSP ) ) ].',...
        tvecInterp_Shim_MPSP_slew_i) ).';

    Shim( :, MPSP_Shim_slew_f_idx_interp ) = ( interp1(...
        round([ ( opt.tEndMPSP - opt.MPSPShimSlewTime ), tvecOrig( Shim_MPSP_slew_f_idxs ), (opt.tEndMPSP) ], sigDigitTime, 'significant' ),...
        [ ( Shim_mag_MPSP .* cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPShimSlewTime ) + Shim_ph_MPSP ) ),...
        -( ( Shim_mag_MPSP .* cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPShimSlewTime ) + Shim_ph_MPSP ) ) / opt.MPSPShimSlewTime ) * ( tvecOrig( Shim_MPSP_slew_f_idxs ) - opt.tEndMPSP ),...
        zeros( opt.numZCoils, 1 ) ].',...
        tvecInterp_Shim_MPSP_slew_f) ).';

    % constant MPSP waveforms
    Shim( :, MPSP_Shim_const_idx_interp ) = ( interp1(...
        round([ ( opt.tStMPSP + opt.MPSPShimSlewTime ), tvecOrig( Shim_MPSP_wave_idxs ), (opt.tEndMPSP - opt.MPSPShimSlewTime) ], sigDigitTime, 'significant' ),...
        (Shim_mag_MPSP .* cos( wz_mpptx * (  [ ( opt.tStMPSP + opt.MPSPShimSlewTime ), tvecOrig( Shim_MPSP_wave_idxs ), (opt.tEndMPSP - opt.MPSPShimSlewTime) ] - opt.tStMPSP ) + Shim_ph_MPSP )).',...
        tvecInterp_MPSP_Shim_const) ).';

else
    Shim = zeros( 1, ( numTimePoints ) );
end

%% Assign values to waveform struct
waveform = struct;

waveform.gyro = 267.5e6;

waveform.RF = RF;
waveform.RFphasor = RFphasor;
waveform.Freq = Freq;
waveform.breal_ORSP_rshp = breal_ORSP_rshp;
waveform.bimag_ORSP_rshp = bimag_ORSP_rshp;
waveform.breal_MPSP_rshp = breal_MPSP_rshp;
waveform.bimag_MPSP_rshp = bimag_MPSP_rshp;

waveform.Grad = Grad;
waveform.Gx_ORSP = Gx_ORSP;
waveform.Gy_ORSP = Gy_ORSP;
waveform.Gz_ORSP = Gz_ORSP;
waveform.G_Blip = G_Blip;
waveform.G_real_MPSP = G_real_MPSP;
waveform.G_imag_MPSP = G_imag_MPSP;
waveform.G_mag_MPSP = G_mag_MPSP;
waveform.G_ph_MPSP = G_ph_MPSP;

waveform.Shim = Shim;
if numZCoils > 0
    waveform.shim_Blip = Shim_Blip;
    waveform.shimreal_MPSP = Shimreal_MPSP;
    waveform.shimimag_MPSP = Shimimag_MPSP;
    waveform.shim_mag_MPSP = Shim_mag_MPSP;
    waveform.shim_ph_MPSP = Shim_ph_MPSP;
end

waveform.numZCoils = numZCoils;
waveform.numXYCoils = numXYCoils;

waveform.numPos = opt.numPos;
waveform.numVars = opt.numVars;

waveform.dwxyvec = zeros( 1, numTimePoints );
waveform.constantRotatingFrame = opt.constantRotatingFrame;

waveform.numTimePoints = numTimePoints;
waveform.tvec = tvecInterp;
waveform.dt = dtInterp;
waveform.dtvec = dtvecInterp;
waveform.pulseLength = pulseLength;

waveform.tStORSP = opt.tStORSP;
waveform.tEndORSP = opt.tEndORSP;
waveform.tStBlip = opt.tStBlip;
waveform.tEndBlip = opt.tEndBlip;
waveform.tStMPSP = opt.tStMPSP;
waveform.tEndMPSP = opt.tEndMPSP;

waveform.BlipGradSlewTime = opt.BlipGradSlewTime;
waveform.BlipShimSlewTime = opt.BlipShimSlewTime;
waveform.MPSPGradSlewTime = opt.MPSPGradSlewTime;
waveform.MPSPShimSlewTime = opt.MPSPShimSlewTime;

waveform.dwxy_mpptx = dwxy_mpptx;
waveform.dfxy_mpptx = dfxy_mpptx;
waveform.wz_mpptx = wz_mpptx;
waveform.fz_mpptx = fz_mpptx;

waveform.Z0 = opt.Z0;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %