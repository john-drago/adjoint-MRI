function [ waveform ] = generateSPINSWaveform_base( p, opt )
% [ RF, Grad, Shim, RFphasor, Freq, tvec, params ] fields of waveform
% Generate waveform for the hard pulse integration.

% Get unscaled parameters -- currently expecting unscaled parameters to be
% passed
% scVec = opt.scVec;
% p = pSc( : ) .* scVec;
p = p( : );

%% Get SPINS parameters
kmax = opt.min_kmax + ( p( opt.kmax_idx ) + opt.scVec( opt.kmax_idx ) );
a = opt.min_a + ( p( opt.a_idx ) + opt.scVec( opt.a_idx ) );
b = opt.min_b + ( p( opt.b_idx ) + opt.scVec( opt.b_idx ) );
u = opt.min_u + ( p( opt.u_idx ) + opt.scVec( opt.u_idx ) );
v = opt.min_v + ( p( opt.v_idx ) + opt.scVec( opt.v_idx ) );

s = struct;
s.kmax = kmax;
s.a = a;
s.b = b;
s.u = u;
s.v = v;
s.T = opt.Tspins;
s.gyro = opt.gyro;

initSlewTime = opt.initSlewTime;
init_slew_i = opt.init_slew_i;
init_slew_f = opt.init_slew_f;
% init_slew_sc = opt.init_slew_sc;
% init_slew_num = opt.init_slew_num;

finalSlewTime = opt.finalSlewTime;
final_slew_i = opt.final_slew_i;
final_slew_f = opt.final_slew_f;
% final_slew_sc = opt.final_slew_sc;
% final_slew_num = opt.final_slew_num;

% spins_i = opt.spins_i;
% spins_f = opt.spins_f;
% spins_int_i = opt.spins_int_i;
% spins_int_f = opt.spins_int_f;
spins_int_idx = opt.spins_int_idx;

%% Initialize Waveforms
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
breal_rshp = transpose( reshape( breal_vec, [ numTimePoints, numXYCoils ] ) );

% bimag
bimag_idx = opt.bimag_idx;
bimag_vec = p( bimag_idx );
bimag_rshp = transpose( reshape( bimag_vec, [ numTimePoints, numXYCoils ] ) );

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
grad_slew_init = ( Gxyz_fn( 0, s ) / initSlewTime ) * opt.tvec( init_slew_i:init_slew_f );
grad_spins_int = Gxyz_fn( opt.tvec( spins_int_idx ) - initSlewTime, s );
grad_slew_final = -( Gxyz_fn( s.T - finalSlewTime, s ) / finalSlewTime ) * ( opt.tvec( final_slew_i:final_slew_f ) - opt.pulseLength );

Grad( :, init_slew_i:init_slew_f ) = grad_slew_init;
Grad( :, spins_int_idx ) = grad_spins_int;
Grad( :, final_slew_i:final_slew_f ) = grad_slew_final;

%% Add Shim points
if opt.numZCoils > 0
    Shim = zeros( numZCoils, numTimePoints );

    % shim
    shim_idx = opt.shim_idx;
    shim_vec = p( shim_idx );
    shim_rshp = transpose( reshape( shim_vec, [ numTimePoints, numZCoils ] ) );

    if strcmpi( opt.structtype, "opt" )

        Shim( :, 1:end ) = shim_rshp;

    elseif strcmpi( opt.structtype, "val" )

        Shim( :, 1:end ) = interp1( [ 0, opt.opt_tvec, opt.pulseLength ].', [ zeros( numZCoils, 1), shim_rshp, zeros( numZCoils, 1) ].', opt.tvec.' ).';
    end
else
    Shim = zeros( 1, numTimePoints );
end

%% Create waveform struct
waveform = generateWaveformBase( breal, bimag, brealphasor, bimagphasor,...
    Grad, Shim, tvec, dtvec, dwxyvec, opt);

waveform.p = p;

waveform.pulseLength = opt.pulseLength;

waveform.kmax = kmax;
waveform.a = a;
waveform.b = b;
waveform.u = u;
waveform.v = v;
waveform.T = s.T;
waveform.s = s;
waveform.Tspins = opt.Tspins;
waveform.initSlewTime = initSlewTime;
waveform.finalSlewTime = finalSlewTime;

end