function opt = determinePercentConstraints( val, opt, pulse )

%% Generate waveforms
wv = opt.generateWaveforms( opt.pOpt, opt );
wv = calculateWaveformSlewAccel( wv );

%% initialize constraint v
constraintvals = dictionary;

%% iterate over pulse constraints
constKeys = keys(pulse.constraints);

for cc = 1:length( constKeys )
    
    currConst = string( constKeys( cc ) );
    currConstVal = pulse.constraints( currConst );

    if matches( currConst, [ "total-RF-power"; "RF-power" ], 'IgnoreCase', true )
        constraintvals( currConst ) = val.totalRFPower / currConstVal;
    elseif matches( currConst, "max-RF-power", 'IgnoreCase', true )
        constraintvals( currConst ) = val.maxRFPower / currConstVal;
    elseif matches( currConst, "RF-slew-rate", 'IgnoreCase', true )
        constraintvals( currConst ) = max(...
            [ val.breal_phasor_slew_max, val.bimag_phasor_slew_max ] ) / currConstVal;
    elseif matches( currConst, "RF-accel", 'IgnoreCase', true )
        constraintvals( currConst ) = max(...
            [ val.breal_phasor_accel_max, val.bimag_phasor_accel_max ] ) / currConstVal;
    elseif matches( currConst, "RF-max", 'IgnoreCase', true )
        constraintvals( currConst ) = max(...
            sqrt( wv.brealphasor.^2 + wv.bimagphasor.^2 ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "RF-bandwidth", 'IgnoreCase', true )
        RFcomp = complex( wv.breal, wv.bimag );
        fftRFcomp = fftshift( fft( RFcomp, [], 2 ) / opt.numTimePoints, 2 ); 
        fvecRFcomp = fvecDFT( 1 / ( wv.tvec(2) - wv.tvec(1) ), length( wv.tvec ), true );
        outOfBandwidth = ( abs( fvecRFcomp ) >= ( val.RFBandwidth_band + 1e-6 ) );
        constraintvals( currConst ) = max(...
            sum( opt.adjdf * abs( fftRFcomp( :, outOfBandwidth ) ), 2 ),...
            [], 'all' ) / opt.RFBandwidth_constr;
    elseif matches( currConst, "peak-local-SAR", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakLocalSAR / currConstVal;
    elseif matches( currConst, "peak-global-SAR", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakGlobalSAR / currConstVal;
    elseif matches( currConst, [ "average-local-SAR"; "avg-local-SAR" ], 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakAvgLocalSAR / currConstVal;
    elseif matches( currConst, [ "average-global-SAR"; "avg-global-SAR" ], 'IgnoreCase', true )
        constraintvals( currConst ) = val.avgGlobalSAR / currConstVal;
     elseif matches( currConst, "average-E-10m", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakAvgE10m / currConstVal;
    elseif matches( currConst, "average-H-10m", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakAvgH10m / currConstVal;
    elseif matches( currConst, "peak-E-10m", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakE10m / currConstVal;
    elseif matches( currConst, "peak-H-10m", 'IgnoreCase', true )
        constraintvals( currConst ) = val.peakH10m / currConstVal;
    elseif matches( currConst, "grad-max", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Grad ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "grad-max-MPSP", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Grad( :, opt.MPSP_i:opt.MPSP_f ) ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "grad-blip-max", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Grad( :, opt.Blip_i:opt.Blip_f ) ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "grad-slew-rate", 'IgnoreCase', true )
        constraintvals( currConst ) = val.Grad_slew_max / currConstVal;
    elseif matches( currConst, "grad-accel", 'IgnoreCase', true )
        constraintvals( currConst ) = val.Grad_accel_max / currConstVal;
    elseif matches( currConst, "shim-max", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Shim ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "shim-max-MPSP", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Shim( :, opt.MPSP_i:opt.MPSP_f ) ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "shim-blip-max", 'IgnoreCase', true )
        constraintvals( currConst ) = max( abs( wv.Shim( :, opt.Blip_i:opt.Blip_f ) ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "shim-total", 'IgnoreCase', true )
        constraintvals( currConst ) = max( sum( abs( wv.Shim ), 1 ), [], 'all' ) / currConstVal;
    elseif matches( currConst, "shim-slew-rate", 'IgnoreCase', true )
        constraintvals( currConst ) = val.Shim_slew_max / currConstVal;
    elseif matches( currConst, "shim-accel", 'IgnoreCase', true )
        constraintvals( currConst ) = val.Shim_accel_max / currConstVal;
    elseif matches( currConst, "max-pulse-length", 'IgnoreCase', true )
        opt = opt.get_tvec( opt.pOpt, opt );
        constraintvals( currConst ) = opt.pulseLength / currConstVal;
    elseif matches( currConst, "min-pulse-length", 'IgnoreCase', true )
        opt = opt.get_tvec( opt.pOpt, opt );
        constraintvals( currConst ) = opt.pulseLength / currConstVal;
    else
        error( "Don't know how to process this constraint" );
    end

end

%% Assign to opt struct
opt.constraintvals = constraintvals;

end