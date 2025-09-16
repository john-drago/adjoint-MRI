function wv = calculateWaveformSlewAccel( wv )

wv.breal_slew = calculateSlew( wv.dtvec, wv.tvec, wv.breal );
wv.breal_accel = calculateAccel( wv.dtvec, wv.tvec, wv.breal );

wv.bimag_slew = calculateSlew( wv.dtvec, wv.tvec, wv.bimag );
wv.bimag_accel = calculateAccel( wv.dtvec, wv.tvec, wv.bimag );

[ wv.brealphasor_slew, wv.dtvec_slew, wv.tvec_slew ] = calculateSlew( wv.dtvec, wv.tvec, wv.brealphasor );
[ wv.brealphasor_accel, wv.dtvec_accel, wv.tvec_accel ] = calculateAccel( wv.dtvec, wv.tvec, wv.brealphasor );

wv.bimagphasor_slew = calculateSlew( wv.dtvec, wv.tvec, wv.bimagphasor );
wv.bimagphasor_accel = calculateAccel( wv.dtvec, wv.tvec, wv.bimagphasor );

wv.Grad_slew = calculateSlew( wv.dtvec, wv.tvec, wv.Grad );
wv.Grad_accel = calculateAccel( wv.dtvec, wv.tvec, wv.Grad );

wv.Shim_slew = calculateSlew( wv.dtvec, wv.tvec, wv.Shim );
wv.Shim_accel = calculateAccel( wv.dtvec, wv.tvec, wv.Shim );

end