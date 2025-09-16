function [ opt, wv ] = gpuArrayAdjointPWC_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

% opt.tvec = gpuArray( opt.tvec );
% opt.dtvec = gpuArray( opt.dtvec );
% opt.gyro = gpuArray( opt.gyro );

end