function partialfpartialshim = calcpartialfpartialshim( multiplyFactor, partialfpartialBz, opt )
partialfpartialshim = partialfpartialBz .*  reshape( multiplyFactor, [ opt.numPos, 1, opt.numZCoils ] );
end