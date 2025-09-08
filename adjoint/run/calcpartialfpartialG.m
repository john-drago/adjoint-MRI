function partialfpartialG = calcpartialfpartialG( multiplyFactor, partialfpartialBz, opt )
partialfpartialG = partialfpartialBz .* reshape( multiplyFactor, [ opt.numPos, 1, 3 ] );
end