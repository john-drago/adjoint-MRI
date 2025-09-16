function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientPWC_base( partialfpartialp, initSt, wv, opt, tt )
% This function will calculate the gradp_fArray at a given time interval.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mtt, Mttp1, tt );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get given time point
% ti = opt.tvec( tt );
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
grad_opt_idxs = tt : opt.numTimePoints : 3*opt.numTimePoints;
grad_idxs = opt.grad_idx( grad_opt_idxs );
varInfo = addVarIdxs( grad_idxs, varInfo );

partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    calcpartialfpartialG( wv.pos, initSt.partialfpartialBz, opt ),...
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