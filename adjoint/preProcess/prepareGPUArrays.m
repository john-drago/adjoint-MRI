function opt = prepareGPUArrays( opt )

opt.gpu_M0 = gpuArray( opt.M0 );

opt.gpu_b1preal = gpuArray( opt.b1preal );
opt.gpu_b1pimag = gpuArray( opt.b1pimag );
opt.gpu_pos = gpuArray( opt.pos );
opt.gpu_bzsens = gpuArray( opt.bzsens );
opt.gpu_db0 = gpuArray( opt.db0 );

opt.gpu_numZCoils = gpuArray( opt.numZCoils );
opt.gpu_numXYCoils = gpuArray( opt.numXYCoils );
opt.gpu_numPos = gpuArray( opt.numPos );
opt.gpu_gyro = gpuArray( opt.gyro );

% opt.Mtarg = gpuArray( opt.Mtarg );

% opt.gpu_dtvec = gpuArray( opt.dtvec );
% opt.gpu_tvec = gpuArray( opt.tvec );

end