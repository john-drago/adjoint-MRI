function [ opt, wv ] = gpuArrayAdjointFourier_cheb( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

wv.numFourier_RF = gpuArray( wv.numFourier_RF );
wv.numCheb_grad = gpuArray( wv.numCheb_grad );
if opt.numZCoils > 0
    wv.numCheb_shim = gpuArray( wv.numCheb_shim );
end

wv.FBRF = gpuArray( wv.FBRF );

wv.Tn = gpuArray( wv.Tn );
% wv.Tn_rep = gpuArray( wv.Tn_rep );
wv.tdom = gpuArray( wv.tdom );


end