function [ partialfpartialGreal, partialfpartialGimag ] = calculatepartialfpartialG_MPSP(...
    partialfpartialBz, partBpartGmagi, partBpartGphi, wv, opt )

% calculate contribution from Gireal_MPSP toward Gmag and Gphi
partGmagipartGreali = wv.G_real_MPSP ./ wv.G_mag_MPSP;
partGphipartGreali = -wv.G_imag_MPSP ./ wv.G_mag_MPSP.^2;

% calculate contribution from Giimag_MPSP toward Gmag and Gphi
partGmagipartGimagi = wv.G_imag_MPSP ./ wv.G_mag_MPSP;
partGphipartGimagi =  wv.G_real_MPSP ./ wv.G_mag_MPSP.^2;

% Catch if there is no magnitude
zeromagIdx = wv.G_mag_MPSP < eps( 1e2 );

if any( zeromagIdx )
    partGmagipartGreali( zeromagIdx ) = 1;
    partGphipartGreali( zeromagIdx ) = 0;
    partGmagipartGimagi( zeromagIdx ) = 0;
    partGphipartGimagi( zeromagIdx ) = 1;
end


% calculate contribution from Gireal_MPSP to B
partialfpartialGreal = calcpartialfpartialG( ...
    ( ( wv.pos ) .* repmat(  ( partBpartGmagi .* partGmagipartGreali + partBpartGphi .* partGphipartGreali ).', [ opt.numPos, 1 ]  ) ),...
    partialfpartialBz, opt );

% calculate contribution from Giimag_MPSP to B
partialfpartialGimag = calcpartialfpartialG( ...
    ( ( wv.pos ) .* repmat( ( partBpartGmagi .* partGmagipartGimagi + partBpartGphi .* partGphipartGimagi ).', [ opt.numPos, 1 ] )  ),...
    partialfpartialBz, opt );

end