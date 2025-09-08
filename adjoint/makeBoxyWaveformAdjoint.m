function [ tvecbxy, wvfrmbxy ] = makeBoxyWaveformAdjoint( tvec, dtvec, wvfrm )
% This function will make a boxy waveform from a wvfrm assumed to be
% sampled in the center of dt periods. The dt lengths are encoded in the
% dtvec, and the tvec tells the where the samplings occur.

tvecbxy = zeros( size( tvec, 1 ), 2 * size( tvec, 2 ) );
tvecbxy( 1, 1:2:end ) = tvec - dtvec/2;
tvecbxy( 1, 2:2:end ) = tvec + dtvec/2;
tvecbxy = [ 0, tvecbxy, tvec(end) + dtvec(end)/2 ];

wvfrmbxy = zeros( size( wvfrm, 1 ), 2 * size( wvfrm, 2 ) );
wvfrmbxy( :, 1:2:end ) = wvfrm;
wvfrmbxy( :, 2:2:end ) = wvfrm;
wvfrmbxy = [ zeros( size( wvfrm, 1 ), 1 ), wvfrmbxy, zeros( size( wvfrm, 1 ), 1 ) ];

end