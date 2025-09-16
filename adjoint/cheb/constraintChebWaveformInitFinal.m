function [ Awv_if, bwv_if ] =  constraintChebWaveformInitFinal(...
    wv_idx, numCheb, opt )

wv_idx_rshp = reshape( wv_idx, numCheb, [ ] );
numWv = size( wv_idx_rshp, 2 );
Awv_if = zeros( 2*numWv, opt.numVars );

Tif = evalChebClenshaw( [ 0; opt.pulseLength ], eye( numCheb ), opt.tdom );
Ti = repmat( Tif( 1, : ), [ numWv, 1 ] );
Tf = repmat( Tif( 2, : ), [ numWv, 1 ] );

Airidx = repmat( transpose(1:numWv), [ 1 numCheb ]);
Afridx = repmat( transpose(numWv+(1:numWv)), [ 1 numCheb ]);
Acidx = transpose( wv_idx_rshp );

linidx = sub2ind( size( Awv_if ), Airidx( : ), Acidx( : ) );
linfdx = sub2ind( size( Awv_if ), Afridx( : ), Acidx( : ) );

Awv_if( linidx ) = Ti( : );
Awv_if( linfdx ) = Tf( : );

bwv_if = zeros( 2*numWv, 1 );

end