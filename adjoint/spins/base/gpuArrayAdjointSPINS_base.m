function [ opt, wv ] = gpuArrayAdjointSPINS_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

% opt.tvec = gpuArray( opt.tvec );
% opt.dtvec = gpuArray( opt.dtvec );
% opt.gyro = gpuArray( opt.gyro );

% wv.pulseLength = gpuArray( wv.pulseLength );
% wv.finalSlewTime = gpuArray( wv.finalSlewTime );
% wv.initSlewTime = gpuArray( wv.initSlewTime );
% wv.kmax = gpuArray( wv.kmax );
% wv.a = gpuArray( wv.a );
% wv.b = gpuArray( wv.b );
% wv.u = gpuArray( wv.u );
% wv.v = gpuArray( wv.v );
% wv.T = gpuArray( wv.T );
% wv.Tspins = gpuArray( wv.Tspins );
% wv.s.kmax = gpuArray( wv.kmax );
% wv.s.a = gpuArray( wv.s.a );
% wv.s.b = gpuArray( wv.s.b );
% wv.s.u = gpuArray( wv.s.u );
% wv.s.v = gpuArray( wv.s.v );
% wv.s.T = gpuArray( wv.s.T );
% wv.s.gyro = gpuArray( wv.s.gyro );

end