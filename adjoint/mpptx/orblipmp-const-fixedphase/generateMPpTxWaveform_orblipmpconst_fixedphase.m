function [ waveform ] = generateMPpTxWaveform_orblipmpconst_fixedphase( p, opt )
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

RFSlewIntPtNum = opt.RFSlewIntPtNum;
rfSlewScale = ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );

brealphasor = zeros( numXYCoils, numTimePoints );
bimagphasor = zeros( numXYCoils, numTimePoints );
breal = zeros( numXYCoils, numTimePoints );
bimag = zeros( numXYCoils, numTimePoints );
Grad = zeros( 3, numTimePoints );

% Get constants
dwxy_mpptx = opt.dwxy_mpptx;
wz_mpptx = opt.wz_mpptx;
dfxy_mpptx = opt.dfxy_mpptx;
fz_mpptx = opt.fz_mpptx;

%% Add phasor points for RF
% on-resonance subpulse
breal_ORSP_p_idx = opt.breal_ORSP_idx;
bimag_ORSP_p_idx = opt.bimag_ORSP_idx;

ORSP_i = opt.ORSP_i;
ORSP_f = opt.ORSP_f;
ORSP_RF_Slew_i = opt.ORSP_RF_Slew_i;
ORSP_RF_Slew_f = opt.ORSP_RF_Slew_f;

ORSP_idx = ORSP_i : ORSP_f;
ORSP_Slew_i_idx = ORSP_i : (ORSP_RF_Slew_i - 1);
ORSP_Slew_f_idx = (ORSP_RF_Slew_f + 1) : ORSP_f;
ORSP_Const_idx = ORSP_RF_Slew_i : ORSP_RF_Slew_f;

breal_ORSP = p( breal_ORSP_p_idx );
bimag_ORSP = p( bimag_ORSP_p_idx );

if strcmpi( opt.structtype, "opt" )

    brealphasor( :, ORSP_Const_idx ) = repmat( breal_ORSP, [ 1, length( ORSP_Const_idx ) ] );
    bimagphasor( :, ORSP_Const_idx ) = repmat( bimag_ORSP, [ 1, length( ORSP_Const_idx ) ] );
    brealphasor( :, ORSP_Slew_i_idx ) = breal_ORSP * rfSlewScale;
    bimagphasor( :, ORSP_Slew_i_idx ) = bimag_ORSP * rfSlewScale;
    brealphasor( :, ORSP_Slew_f_idx ) = fliplr( brealphasor( :, ORSP_Slew_i_idx ) );
    bimagphasor( :, ORSP_Slew_f_idx ) = fliplr( bimagphasor( :, ORSP_Slew_i_idx ) );

    breal( :, ORSP_idx ) = brealphasor( :, ORSP_idx );
    bimag( :, ORSP_idx ) = bimagphasor( :, ORSP_idx );

elseif strcmpi( opt.structtype, "val" )

    brealphasor( :, ORSP_Slew_i_idx ) = ( breal_ORSP / opt.RFSlewTime ) *...
        ( tvec( ORSP_Slew_i_idx ) - opt.tStORSP );
    brealphasor( :, ORSP_Const_idx ) = repmat( breal_ORSP, [ 1, length( ORSP_Const_idx ) ] );
    brealphasor( :, ORSP_Slew_f_idx ) = -( breal_ORSP / opt.RFSlewTime ) *...
        ( tvec( ORSP_Slew_f_idx ) - opt.tEndORSP );

    bimagphasor( :, ORSP_Slew_i_idx ) = ( bimag_ORSP / opt.RFSlewTime ) *...
        ( tvec( ORSP_Slew_i_idx ) - opt.tStORSP );
    bimagphasor( :, ORSP_Const_idx ) = repmat( bimag_ORSP, [ 1, length( ORSP_Const_idx ) ] );
    bimagphasor( :, ORSP_Slew_f_idx ) = -( bimag_ORSP / opt.RFSlewTime ) *...
        ( tvec( ORSP_Slew_f_idx ) - opt.tEndORSP );

    breal( :, ORSP_idx ) = brealphasor( :, ORSP_idx );
    bimag( :, ORSP_idx ) = bimagphasor( :, ORSP_idx );

end

% multiphoton subpulse
breal_MPSP_p_idx = opt.breal_MPSP_idx;
bimag_MPSP_p_idx = opt.bimag_MPSP_idx;

MPSP_i = opt.MPSP_i;
MPSP_f = opt.MPSP_f;
MPSP_RF_Slew_i = opt.MPSP_RF_Slew_i;
MPSP_RF_Slew_f = opt.MPSP_RF_Slew_f;

MPSP_idx = MPSP_i : MPSP_f;
MPSP_Slew_i_idx = MPSP_i : (MPSP_RF_Slew_i - 1);
MPSP_Slew_f_idx = (MPSP_RF_Slew_f + 1) : MPSP_f;
MPSP_Const_idx = MPSP_RF_Slew_i : MPSP_RF_Slew_f;

breal_MPSP = p( breal_MPSP_p_idx );
bimag_MPSP = p( bimag_MPSP_p_idx );

cosdwxy = cos( dwxy_mpptx * ( tvec( MPSP_idx ) - opt.tStMPSP) );
sindwxy = sin( dwxy_mpptx * ( tvec( MPSP_idx ) - opt.tStMPSP) );

if strcmpi( opt.structtype, "opt" )

    brealphasor( :, MPSP_Const_idx ) = repmat( breal_MPSP, [ 1, length( MPSP_Const_idx ) ] );
    bimagphasor( :, MPSP_Const_idx ) = repmat( bimag_MPSP, [ 1, length( MPSP_Const_idx ) ] );
    brealphasor( :, MPSP_Slew_i_idx ) = breal_MPSP * rfSlewScale;
    bimagphasor( :, MPSP_Slew_i_idx ) = bimag_MPSP * rfSlewScale;
    brealphasor( :, MPSP_Slew_f_idx ) = fliplr( brealphasor( :, MPSP_Slew_i_idx ) );
    bimagphasor( :, MPSP_Slew_f_idx ) = fliplr( bimagphasor( :, MPSP_Slew_i_idx ) );

    breal( :, MPSP_idx ) = ( brealphasor( :, MPSP_idx ) .* cosdwxy - bimagphasor( :, MPSP_idx ) .* sindwxy );
    bimag( :, MPSP_idx ) = ( bimagphasor( :, MPSP_idx ) .* cosdwxy + brealphasor( :, MPSP_idx ) .* sindwxy );

elseif strcmpi( opt.structtype, "val" )

    brealphasor( :, MPSP_Slew_i_idx ) = ( breal_MPSP / opt.RFSlewTime ) *...
        ( tvec( MPSP_Slew_i_idx ) - opt.tStMPSP );
    brealphasor( :, MPSP_Const_idx ) = repmat( breal_MPSP, [ 1, length( MPSP_Const_idx ) ] );
    brealphasor( :, MPSP_Slew_f_idx ) = -( breal_MPSP / opt.RFSlewTime ) *...
        ( tvec( MPSP_Slew_f_idx ) - opt.tEndMPSP );

    bimagphasor( :, MPSP_Slew_i_idx ) = ( bimag_MPSP / opt.RFSlewTime ) *...
        ( tvec( MPSP_Slew_i_idx ) - opt.tStMPSP );
    bimagphasor( :, MPSP_Const_idx ) = repmat( bimag_MPSP, [ 1, length( MPSP_Const_idx ) ] );
    bimagphasor( :, MPSP_Slew_f_idx ) = -( bimag_MPSP / opt.RFSlewTime ) *...
        ( tvec( MPSP_Slew_f_idx ) - opt.tEndMPSP );

    breal( :, MPSP_idxs ) = ( brealphasor( :, MPSP_idx ) .* cosdwxy - bimagphasor( :, MPSP_idx ) .* sindwxy );
    bimag( :, MPSP_idxs ) = ( bimagphasor( :, MPSP_idx ) .* cosdwxy + brealphasor( :, MPSP_idx ) .* sindwxy );

end

%% Add Grad points
% blip period
Gx_Blip_p_idx = opt.Gx_Blip_idx;
Gy_Blip_p_idx = opt.Gy_Blip_idx;
Gz_Blip_p_idx = opt.Gz_Blip_idx;

Gx_Blip = p( Gx_Blip_p_idx );
Gy_Blip = p( Gy_Blip_p_idx );
Gz_Blip = p( Gz_Blip_p_idx );

Grad_Blip_slew_i_idxs = opt.Blip_i : ( opt.Blip_Grad_Slew_i-1 );
Grad_Blip_slew_f_idxs = ( opt.Blip_Grad_Slew_f+1 ) : opt.Blip_f;
Grad_Blip_const_idxs = opt.Blip_Grad_Slew_i : opt.Blip_Grad_Slew_f;

G_Blip = [ Gx_Blip; Gy_Blip; Gz_Blip ];

Grad( :, Grad_Blip_slew_i_idxs ) = ( G_Blip / opt.BlipGradSlewTime ) *...
    ( tvec( Grad_Blip_slew_i_idxs ) - opt.tStBlip );
Grad( :, Grad_Blip_slew_f_idxs ) = -( G_Blip / opt.BlipGradSlewTime ) *...
    ( tvec( Grad_Blip_slew_f_idxs ) - opt.tEndBlip );

Grad( :, Grad_Blip_const_idxs ) = repmat( G_Blip, [ 1 length( Grad_Blip_const_idxs ) ] );

% multiphoton subpulse
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

Grad( :, G_MPSP_slew_i_idxs ) = ( ( G_mag_MPSP / opt.MPSPGradSlewTime ) .*...
    cos( wz_mpptx * opt.MPSPGradSlewTime +  G_ph_MPSP ) ) *...
    ( tvec( G_MPSP_slew_i_idxs ) - opt.tStMPSP );
Grad( :, G_MPSP_slew_f_idxs ) = -( ( G_mag_MPSP / opt.MPSPGradSlewTime ) .*...
    cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPGradSlewTime ) +  G_ph_MPSP ) ) *...
    ( tvec( G_MPSP_slew_f_idxs ) - opt.tEndMPSP );
Grad( :, G_MPSP_wave_idxs ) = G_mag_MPSP .* ( cos( wz_mpptx * ( tvec( G_MPSP_wave_idxs ) - opt.tStMPSP ) + G_ph_MPSP ) );

%% Add Shim points
if numZCoils > 0
    Shim = zeros( numZCoils, numTimePoints );

    % multiphoton subpulse
    shimmag_MPSP_p_idx = opt.shimmag_MPSP_idx;
    shimph_MPSP_p_idx = opt.shimph_MPSP_idx;

    shim_mag_MPSP = p( shimmag_MPSP_p_idx );
    shim_ph_MPSP = p( shimph_MPSP_p_idx );

    shim_MPSP_slew_i_idx = opt.MPSP_i:( opt.MPSP_Shim_Slew_i-1 );
    shim_MPSP_slew_f_idx = ( opt.MPSP_Shim_Slew_f+1 ):( opt.MPSP_f );
    shim_MPSP_wave_idx = ( opt.MPSP_Shim_Slew_i ) : opt.MPSP_Shim_Slew_f;

    % slew waveforms
    Shim( :, shim_MPSP_slew_i_idx ) = ( ( shim_mag_MPSP / opt.MPSPShimSlewTime ) .*...
        cos( wz_mpptx * opt.MPSPShimSlewTime +  shim_ph_MPSP ) )*...
        ( tvec( shim_MPSP_slew_i_idx ) - opt.tStMPSP );
    Shim( :, shim_MPSP_slew_f_idx ) = -( ( shim_mag_MPSP / opt.MPSPShimSlewTime ) .*...
        cos( wz_mpptx * ( opt.tEndMPSP - opt.tStMPSP - opt.MPSPShimSlewTime ) +  shim_ph_MPSP ) )*...
        ( tvec( shim_MPSP_slew_f_idx ) - opt.tEndMPSP );
    Shim( :, shim_MPSP_wave_idx ) = shim_mag_MPSP .* ( cos( wz_mpptx * ( tvec( shim_MPSP_wave_idx ) - opt.tStMPSP ) + shim_ph_MPSP ) );

    % blip period
    shim_Blip_scale_idx = opt.shimblipscale_idx;

    shim_Blip_scale = p( shim_Blip_scale_idx );

    shim_Blip_slew_i_idxs = opt.Blip_i:( opt.Blip_Shim_Slew_i-1 );
    shim_Blip_slew_f_idxs = ( opt.Blip_Shim_Slew_f+1 ):( opt.Blip_f );
    shim_Blip_const_idxs = ( opt.Blip_Shim_Slew_i ) : opt.Blip_Shim_Slew_f;

    Shim( :, shim_Blip_slew_i_idxs ) = ( ( shim_Blip_scale * shim_mag_MPSP ) / opt.BlipShimSlewTime ) *...
        ( tvec( shim_Blip_slew_i_idxs ) - opt.tStBlip );
    Shim( :, shim_Blip_slew_f_idxs ) = -( ( shim_Blip_scale * shim_mag_MPSP ) / opt.BlipShimSlewTime ) *...
        ( tvec( shim_Blip_slew_f_idxs ) - opt.tEndBlip );

    Shim( :, shim_Blip_const_idxs ) = repmat( ( shim_Blip_scale * shim_mag_MPSP ), [ 1 length( shim_Blip_const_idxs ) ] );
else
    Shim = zeros( 1, numTimePoints );
end

%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.brealphasor = brealphasor;
waveform.bimagphasor = bimagphasor;
waveform.breal = breal;
waveform.bimag = bimag;
waveform.breal_ORSP = breal_ORSP;
waveform.bimag_ORSP = bimag_ORSP;
waveform.breal_MPSP = breal_MPSP;
waveform.bimag_MPSP = bimag_MPSP;

waveform.G_Blip = G_Blip;
waveform.G_real_MPSP = G_real_MPSP;
waveform.G_imag_MPSP = G_imag_MPSP;
waveform.G_mag_MPSP = G_mag_MPSP;
waveform.G_ph_MPSP = G_ph_MPSP;

if numZCoils > 0
    waveform.shim_Blip_scale = shim_Blip_scale;
    waveform.shim_Blip = shim_Blip_scale * shim_mag_MPSP;
    waveform.shim_mag_MPSP = shim_mag_MPSP;
    waveform.shim_ph_MPSP = shim_ph_MPSP;
end

waveform.tORSPIntPtNum = opt.tORSPIntPtNum;
waveform.tMPSPIntPtNum = opt.tMPSPIntPtNum;
waveform.num_ORSP = opt.num_ORSP;
waveform.num_MPSP = opt.num_MPSP;

waveform.RFSlewTime = opt.RFSlewTime;

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
waveform.wz_mpptx = wz_mpptx;
waveform.dfxy_mpptx = dfxy_mpptx;
waveform.fz_mpptx = fz_mpptx;

end