function [ val ] = processAdjointPower( val, wv )
% This function will process the power of an optimized adjoint pulse and
% report relevant metrics.

%% Make RF phasor
wv.RF = complex( wv.breal, wv.bimag );
wv.RFphasor = complex( wv.brealphasor, wv.bimagphasor );

%% Pulse Power
[ val.totalRFPower, val.maxRFPower ] = calcTotalMaxRFPower( ...
    wv.RFphasor, wv.dtvec, val.pulseLength, val.Z0, val.dutyCycle );

%% Local SAR
if isfield( val, 'VOPs' )
    [ val.peakLocalSAR, val.peakAvgLocalSAR ] = calcLocalSAR( ...
        val.VOPs, wv.RFphasor, wv.dtvec, val.pulseLength, val.dutyCycle );
end

%% Global SAR
if isfield( val, 'QGlobal' )
    [ val.peakGlobalSAR, val.avgGlobalSAR ] = calcGlobalSAR( ...
        val.QGlobal, wv.RFphasor, wv.dtvec, val.pulseLength, val.dutyCycle );
end

%% Local SED

%% Global SED

%% E 10 m
if isfield( val, 'VOPs_E10m' )
    [ val.peakE10m, val.peakAvgE10m ] = calcEH10m( ...
        val.VOPs_E10m, wv.RFphasor, wv.dtvec, val.pulseLength, val.dutyCycle );
end

%% H 10 m
if isfield( val, 'VOPs_H10m' )
    [ val.peakH10m, val.peakAvgH10m ] = calcEH10m( ...
        val.VOPs_H10m, wv.RFphasor, wv.dtvec, val.pulseLength, val.dutyCycle );
end

end

