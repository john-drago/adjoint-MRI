function partialfpartialSPINS = calcpartialfpartialSPINSparams( dGxyzdparam, pos, partialfpartialBz, opt )
dBzdparam = ( pos * dGxyzdparam );
partialfpartialSPINS = partialfpartialBz .* reshape( dBzdparam, [ opt.numPos, 1, length( opt.param_idx ) ] );
end