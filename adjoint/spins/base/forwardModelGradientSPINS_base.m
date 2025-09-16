function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientSPINS_base( partialfpartialp, initSt, wv, opt, tt )
% This function will calculate the gradp_fArray at a given time interval.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mtt, Mttp1, tt );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get partialB partial parameter for the given time point 
% determine what phase of the pulse we are in

% get time
ti = opt.tvec( tt );
% dti = opt.dtvec( tt );

%% RF
% Calculate contribution from breal
rf_opt_idxs = tt : opt.numTimePoints : (opt.numXYCoils)*opt.numTimePoints;
breal_idxs = ( opt.breal_idx( rf_opt_idxs )  );
varInfo = addVarIdxs( breal_idxs, varInfo );

partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    initSt.partialfpartialBx .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) +...
    initSt.partialfpartialBy .* reshape( wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] ),...
    varInfo );

% Calculate contribution from bimag
bimag_idxs = ( opt.bimag_idx( rf_opt_idxs ) );
varInfo = addVarIdxs( bimag_idxs, varInfo );

partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    initSt.partialfpartialBx .* reshape( ( - wv.b1pimag ), [ opt.numPos, 1, opt.numXYCoils ] ) +...
    initSt.partialfpartialBy .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ),...
    varInfo );

%% Grad
if ( tt >= opt.spins_int_i ) && ( tt <= opt.spins_int_f )

    dGxyzdparam = dGxyzdparam_fn( ti - opt.initSlewTime, wv.s );

elseif ( tt > opt.spins_int_f )

    slewScale = -( ti - opt.pulseLength ) / opt.finalSlewTime;
    dGxyzdparam = slewScale * dGxyzdparam_fn( opt.Tspins - opt.finalSlewTime, wv.s );

elseif ( tt < opt.spins_int_i )

    slewScale = ( ti ) / wv.initSlewTime;
    dGxyzdparam = slewScale * dGxyzdparam_fn( 0, wv.s );

end

dGxyzdparam = squeeze( dGxyzdparam );

if opt.useGPU
    dGxyzdparam = gpuArray( dGxyzdparam );
end

% Calculate contribution from gradient parameters
varInfo = addVarIdxs( opt.param_idx, varInfo );

partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    calcpartialfpartialSPINSparams( dGxyzdparam, wv.pos, initSt.partialfpartialBz, opt ),...
    varInfo );

%% Shim
if opt.numZCoils > 0
    shim_opt_idxs = tt : opt.numTimePoints : opt.numZCoils * opt.numTimePoints;
    shim_idxs = opt.shim_idx( shim_opt_idxs );
    varInfo = addVarIdxs( shim_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialfpartialshim( wv.bzsens, initSt.partialfpartialBz, opt ),...
        varInfo );
end

end