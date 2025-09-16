function [ partialfpartialGconst_real, partialfpartialGconst_imag ] = calculatepartialfpartialG_MPSP_const(...
    partialfpartialBz, wv, opt, ti )

partBpartGmagi = cos( wv.wz_mpptx * ( ti - wv.tStMPSP ) + wv.G_ph_MPSP );
partBpartGphi = ( - wv.G_mag_MPSP ) .* sin( wv.wz_mpptx * ( ti - wv.tStMPSP ) + wv.G_ph_MPSP );

[ partialfpartialGconst_real, partialfpartialGconst_imag ] = calculatepartialfpartialG_MPSP(...
    partialfpartialBz, partBpartGmagi, partBpartGphi, wv, opt );

end
