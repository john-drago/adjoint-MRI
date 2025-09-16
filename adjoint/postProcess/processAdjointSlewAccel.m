function [ val, wv ] = processAdjointSlewAccel( val, wv )
% This function will process the slew rate and acceleration of the RF,
% Grad, and Shim for an optimized pulse.

%% Initialize slew and accel matrices
wv = calculateWaveformSlewAccel( wv );

%% Find max and min
dim = 'all';
val.breal_phasor_slew_max = max( abs( wv.brealphasor_slew ), [], dim );
val.breal_phasor_accel_max = max( abs( wv.brealphasor_accel ), [], dim );

val.bimag_phasor_slew_max = max( abs( wv.bimagphasor_slew ), [], dim );
val.bimag_phasor_accel_max = max( abs( wv.bimagphasor_accel ), [], dim );

val.Grad_slew_max = max( abs( wv.Grad_slew ), [], dim );
val.Grad_accel_max = max( abs( wv.Grad_accel ), [], dim );

val.Shim_slew_max = max( abs( wv.Shim_slew ), [], dim );
val.Shim_accel_max = max( abs( wv.Shim_accel ), [], dim );

end
