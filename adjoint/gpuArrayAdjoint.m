function wv = gpuArrayAdjoint( wv, opt )

wv.b1preal = gpuArray( wv.b1preal );
wv.b1pimag = gpuArray( wv.b1pimag );
wv.pos = gpuArray( wv.pos );
wv.bzsens = gpuArray( wv.bzsens );
wv.db0 = gpuArray( wv.db0 );

wv.breal = gpuArray( wv.breal );
wv.bimag = gpuArray( wv.bimag );
wv.Grad = gpuArray( wv.Grad );
wv.Shim = gpuArray( wv.Shim );

wv.numZCoils = gpuArray( opt.numZCoils );
wv.numXYCoils = gpuArray( opt.numXYCoils );
wv.numPos = gpuArray( opt.numPos );
wv.numVars = gpuArray( opt.numVars );
wv.gyro = gpuArray( opt.gyro );

wv.dtvec = gpuArray( wv.dtvec );
wv.tvec = gpuArray( wv.tvec );
wv.dwxyvec = gpuArray( wv.dwxyvec ); 

end