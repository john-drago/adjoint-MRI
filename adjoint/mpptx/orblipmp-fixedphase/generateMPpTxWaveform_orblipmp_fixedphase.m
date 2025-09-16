function [ waveform ] = generateMPpTxWaveform_orblipmp_fixedphase( p, opt )
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

ORSP_idxs = opt.ORSP_i : opt.ORSP_f;

breal_ORSP_vec = p( breal_ORSP_p_idx );
bimag_ORSP_vec = p( bimag_ORSP_p_idx );

breal_ORSP_rshp = reshape( breal_ORSP_vec, [ length( ORSP_idxs ), numXYCoils ] ).';
bimag_ORSP_rshp = reshape( bimag_ORSP_vec, [ length( ORSP_idxs ), numXYCoils ] ).';

if strcmpi( opt.structtype, "opt" )

    brealphasor( :, ORSP_idxs ) = breal_ORSP_rshp;
    bimagphasor( :, ORSP_idxs ) = bimag_ORSP_rshp;

    breal( :, ORSP_idxs ) = breal_ORSP_rshp;
    bimag( :, ORSP_idxs ) = bimag_ORSP_rshp;

elseif strcmpi( opt.structtype, "val" )

    brealphasor( :, ORSP_idxs ) = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ zeros( numXYCoils, 1), breal_ORSP_rshp, zeros( numXYCoils, 1) ].', tvec( ORSP_idxs ) ).';
    bimagphasor( :, ORSP_idxs ) = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ zeros( numXYCoils, 1), bimag_ORSP_rshp, zeros( numXYCoils, 1) ].', tvec( ORSP_idxs ) ).';

    breal( :, ORSP_idxs ) = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ zeros( numXYCoils, 1), breal_ORSP_rshp, zeros( numXYCoils, 1) ].', tvec( ORSP_idxs ) ).';
    bimag( :, ORSP_idxs ) = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ zeros( numXYCoils, 1), bimag_ORSP_rshp, zeros( numXYCoils, 1) ].', tvec( ORSP_idxs ) ).';

end

% multiphoton subpulse
breal_MPSP_p_idx = opt.breal_MPSP_idx;
bimag_MPSP_p_idx = opt.bimag_MPSP_idx;

MPSP_idxs = opt.MPSP_i : opt.MPSP_f;

breal_MPSP_vec = p( breal_MPSP_p_idx );
bimag_MPSP_vec = p( bimag_MPSP_p_idx );

breal_MPSP_rshp = reshape( breal_MPSP_vec, [ length( MPSP_idxs ), numXYCoils ] ).';
bimag_MPSP_rshp = reshape( bimag_MPSP_vec, [ length( MPSP_idxs ), numXYCoils ] ).';

cosdwxy = cos( dwxy_mpptx * ( tvec( MPSP_idxs ) - opt.tStMPSP) );
sindwxy = sin( dwxy_mpptx * ( tvec( MPSP_idxs ) - opt.tStMPSP) );

if strcmpi( opt.structtype, "opt" )

    brealphasor( :, MPSP_idxs ) = breal_MPSP_rshp;
    bimagphasor( :, MPSP_idxs ) = bimag_MPSP_rshp;

    breal( :, MPSP_idxs ) = ( breal_MPSP_rshp .* cosdwxy - bimag_MPSP_rshp .* sindwxy );
    bimag( :, MPSP_idxs ) = ( bimag_MPSP_rshp .* cosdwxy + breal_MPSP_rshp .* sindwxy );

elseif strcmpi( opt.structtype, "val" )

    brealphasor( :, MPSP_idxs ) = interp1( [ opt.tStMPSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStMPSP ) & (opt.opt_tvec <= opt.tEndMPSP ) ), opt.tEndMPSP ], [ zeros( numXYCoils, 1), breal_MPSP_rshp, zeros( numXYCoils, 1) ].', tvec( MPSP_idxs ) ).';
    bimagphasor( :, MPSP_idxs ) = interp1( [ opt.tStMPSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStMPSP ) & (opt.opt_tvec <= opt.tEndMPSP ) ), opt.tEndMPSP ], [ zeros( numXYCoils, 1), bimag_MPSP_rshp, zeros( numXYCoils, 1) ].', tvec( MPSP_idxs ) ).';

    breal( :, MPSP_idxs ) = ( brealphasor( :, MPSP_idxs ) .* cosdwxy - bimagphasor( :, MPSP_idxs ) .* sindwxy );
    bimag( :, MPSP_idxs ) = ( bimagphasor( :, MPSP_idxs ) .* cosdwxy + brealphasor( :, MPSP_idxs ) .* sindwxy );

end

%% Add Grad points
% on-resonance subpulse
Gx_ORSP_p_idx = opt.Gx_ORSP_idx;
Gy_ORSP_p_idx = opt.Gy_ORSP_idx;
Gz_ORSP_p_idx = opt.Gz_ORSP_idx;

if strcmpi( opt.structtype, "opt" )

    Gx_ORSP = p( Gx_ORSP_p_idx ).';
    Gy_ORSP = p( Gy_ORSP_p_idx ).';
    Gz_ORSP = p( Gz_ORSP_p_idx ).';

elseif strcmpi( opt.structtype, "val" )

    Gx_ORSP = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ 0, p( Gx_ORSP_p_idx ).', 0 ], tvec( ORSP_idxs ) );
    Gy_ORSP = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ 0, p( Gy_ORSP_p_idx ).', 0 ], tvec( ORSP_idxs ) );
    Gz_ORSP = interp1( [ opt.tStORSP, opt.opt_tvec( (opt.opt_tvec >= opt.tStORSP ) & (opt.opt_tvec <= opt.tEndORSP ) ), opt.tEndORSP ], [ 0, p( Gz_ORSP_p_idx ).', 0 ], tvec( ORSP_idxs ) );

end

Grad( :, ORSP_idxs ) = [ Gx_ORSP; Gy_ORSP; Gz_ORSP ];

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

waveform.breal_ORSP_rshp = breal_ORSP_rshp;
waveform.bimag_ORSP_rshp = bimag_ORSP_rshp;
waveform.breal_MPSP_rshp = breal_MPSP_rshp;
waveform.bimag_MPSP_rshp = bimag_MPSP_rshp;

waveform.Gx_ORSP = Gx_ORSP;
waveform.Gy_ORSP = Gy_ORSP;
waveform.Gz_ORSP = Gz_ORSP;
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

waveform.num_ORSP = opt.num_ORSP;
waveform.num_MPSP = opt.num_MPSP;

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