function [ opt, wv ] = gpuArrayAdjointMPpTx_orblipmpconst( opt, wv )
% This function will make gpuArrays out of parameters that are needed for
% calculation of forward model and gradients for mpptx_orblipmp.

wv = gpuArrayAdjoint( wv, opt );

if opt.numZCoils > 0
    wv.shim_mag_MPSP = gpuArray( wv.shim_mag_MPSP );
    wv.shim_ph_MPSP = gpuArray( wv.shim_ph_MPSP );
    wv.shimreal_MPSP = gpuArray( wv.shimreal_MPSP );
    wv.shimimag_MPSP = gpuArray( wv.shimimag_MPSP );
end
wv.G_mag_MPSP = gpuArray( wv.G_mag_MPSP );
wv.G_ph_MPSP = gpuArray( wv.G_ph_MPSP );
wv.G_real_MPSP = gpuArray( wv.G_real_MPSP );
wv.G_imag_MPSP = gpuArray( wv.G_imag_MPSP );

wv.num_ORSP = gpuArray( wv.num_ORSP );
wv.num_MPSP = gpuArray( wv.num_MPSP );

wv.tStORSP = gpuArray( wv.tStORSP );
wv.tEndORSP = gpuArray( wv.tEndORSP );
wv.tStBlip = gpuArray( wv.tStBlip );
wv.tEndBlip = gpuArray( wv.tEndBlip );
wv.tStMPSP = gpuArray( wv.tStMPSP );
wv.tEndMPSP = gpuArray( wv.tEndMPSP );
wv.RFSlewTime = gpuArray( wv.RFSlewTime );

wv.BlipGradSlewTime = gpuArray( wv.BlipGradSlewTime );
wv.BlipShimSlewTime = gpuArray( wv.BlipShimSlewTime );
wv.MPSPGradSlewTime = gpuArray( wv.MPSPGradSlewTime );
wv.MPSPShimSlewTime = gpuArray( wv.MPSPShimSlewTime );

wv.dwxy_mpptx = gpuArray( wv.dwxy_mpptx );
wv.dfxy_mpptx = gpuArray( wv.dfxy_mpptx );
wv.wz_mpptx = gpuArray( wv.wz_mpptx );
wv.fz_mpptx = gpuArray( wv.fz_mpptx );


end