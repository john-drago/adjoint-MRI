function [ opt, wv ] = gpuArrayAdjointvarkTP_base( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients.

wv = gpuArrayAdjoint( wv, opt );

wv.breal_ktp_rshp = gpuArray( wv.breal_ktp_rshp );
wv.bimag_ktp_rshp = gpuArray( wv.bimag_ktp_rshp );
wv.grad_ktp_rshp = gpuArray( wv.grad_ktp_rshp );
if opt.numZCoils > 0
    wv.shim_ktp_rshp = gpuArray( wv.shim_ktp_rshp );
end

wv.num_kTP = gpuArray( wv.num_kTP );

wv.RFTiming = gpuArray( wv.RFTiming );
wv.RFLength = gpuArray( wv.RFLength );
wv.RFSlewTime = gpuArray( wv.RFSlewTime );
wv.blipLength = gpuArray( wv.blipLength );
wv.blipTiming = gpuArray( wv.blipTiming );
wv.blipCenters = gpuArray( wv.blipCenters );

wv.gradSlewTime = gpuArray( wv.gradSlewTime );
if opt.numZCoils > 0
    wv.shimSlewTime = gpuArray( wv.shimSlewTime );
end

% opt.tvec = gpuArray( opt.tvec );
% opt.dtvec = gpuArray( opt.dtvec );
% opt.gyro = gpuArray( opt.gyro );

end