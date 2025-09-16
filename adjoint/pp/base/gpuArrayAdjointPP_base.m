function [ opt, wv ] = gpuArrayAdjointPP_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

wv.numVarsPerChannel_RF = gpuArray( wv.numVarsPerChannel_RF );
wv.numPPShape_RF = gpuArray( wv.numPPShape_RF );
wv.PPIdxs_RF = gpuArray( wv.PPIdxs_RF );
wv.PPStartStop_RF = gpuArray( wv.PPStartStop_RF );
wv.PPVarIdxs_RF = gpuArray( wv.PPVarIdxs_RF );
wv.varsToTimepoints_RF = gpuArray( wv.varsToTimepoints_RF );
% wv.shapeFnValsTimePoints_RF_rep = gpuArray( wv.shapeFnValsTimePoints_RF_rep );

wv.numVarsPerChannel_grad = gpuArray( wv.numVarsPerChannel_grad );
wv.numPPShape_grad = gpuArray( wv.numPPShape_grad );
wv.PPIdxs_grad = gpuArray( wv.PPIdxs_grad );
wv.PPStartStop_grad = gpuArray( wv.PPStartStop_grad );
wv.PPVarIdxs_grad = gpuArray( wv.PPVarIdxs_grad );
wv.varsToTimepoints_grad = gpuArray( wv.varsToTimepoints_grad );
% wv.shapeFnValsTimePoints_grad_rep = gpuArray( wv.shapeFnValsTimePoints_grad_rep );

if opt.numZCoils > 0
    wv.numVarsPerChannel_shim = gpuArray( wv.numVarsPerChannel_shim );
    wv.numPPShape_shim = gpuArray( wv.numPPShape_shim );
    wv.PPIdxs_shim = gpuArray( wv.PPIdxs_shim );
    wv.PPStartStop_shim = gpuArray( wv.PPStartStop_shim );
    wv.PPVarIdxs_shim = gpuArray( wv.PPVarIdxs_shim );
    wv.varsToTimepoints_shim = gpuArray( wv.varsToTimepoints_shim );
    % wv.shapeFnValsTimePoints_shim_rep = gpuArray( wv.shapeFnValsTimePoints_shim_rep );
end

% opt.tvec = gpuArray( opt.tvec );
% opt.dtvec = gpuArray( opt.dtvec );
% opt.gyro = gpuArray( opt.gyro );
% opt.Tn = gpuArray( opt.Tn );
% opt.tdom = gpuArray( opt.tdom );

end