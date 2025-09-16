function [ opt, wv ] = gpuArrayAdjointCheb_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

wv.numCheb_RF = gpuArray( wv.numCheb_RF );
wv.breal_rshp = gpuArray( wv.breal_rshp );
wv.bimag_rshp = gpuArray( wv.bimag_rshp );

wv.numCheb_grad = gpuArray( wv.numCheb_grad );
wv.grad_rshp = gpuArray( wv.grad_rshp );

if opt.numZCoils > 0
    wv.numCheb_shim = gpuArray( wv.numCheb_shim );
    wv.shim_rshp = gpuArray( wv.shim_rshp );
end

wv.Tn = gpuArray( wv.Tn );
% wv.Tn_rep = gpuArray( wv.Tn_rep );
wv.tdom = gpuArray( wv.tdom );

% opt.tvec = gpuArray( opt.tvec );
% opt.dtvec = gpuArray( opt.dtvec );
% opt.gyro = gpuArray( opt.gyro );
% opt.Tn = gpuArray( opt.Tn );
% opt.tdom = gpuArray( opt.tdom );

end