function [ opt, wv ] = gpuArrayAdjointFourier_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

wv.numFourier_RF = gpuArray( wv.numFourier_RF );
wv.numFourier_grad = gpuArray( wv.numFourier_grad );
if opt.numZCoils > 0
    wv.numFourier_shim = gpuArray( wv.numFourier_shim );
end

wv.FBRF = gpuArray( wv.FBRF );
wv.FBgrad = gpuArray( wv.FBgrad );
if opt.numZCoils > 0
    wv.FBshim = gpuArray( wv.FBshim );
end


end