function partialBpartialGradShape = calcpartialBpartialGradPP( partialfpartialBz, ShapeFnValrep, opt, wv  )

% partial w.r.t grad
numActiveSF = size(ShapeFnValrep,2);

partialBzpartialgradshape = reshape(...
    ShapeFnValrep .* reshape( wv.pos, [ opt.numPos, 1, 3 ] ),...
    [ opt.numPos, 1, numActiveSF * 3 ] );

partialBpartialGradShape =...
    partialfpartialBz .* partialBzpartialgradshape;

% clear partialBzpartialgradcheb;

end