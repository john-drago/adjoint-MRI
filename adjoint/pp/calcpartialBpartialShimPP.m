function partialBpartialShimShape = calcpartialBpartialShimPP( partialfpartialBz, ShapeFnValrep, opt, wv  )

% partial w.r.t shim
numActiveSF = size(ShapeFnValrep,2);

partialBzpartialshimshape = reshape(...
    ShapeFnValrep .* reshape( wv.bzsens, [opt.numPos, 1, opt.numZCoils ] ),...
    [ opt.numPos, 1, numActiveSF*opt.numZCoils ] );

partialBpartialShimShape =...
    partialfpartialBz .* partialBzpartialshimshape;

% clear partialBzpartialshimcheb;

end