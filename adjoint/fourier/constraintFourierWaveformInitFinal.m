function [ Awv_if, bwv_if ] =  constraintFourierWaveformInitFinal(...
    wv_idx_rshp, numFourier, opt )

numWv = size( wv_idx_rshp, 2 );
Awv_if = zeros( numWv, opt.numVars );

Fif = ones( numWv, numFourier );

Aifridx = repmat( transpose(1:numWv), [ 1, numFourier ]);
Acidx = transpose( wv_idx_rshp );

linifdx = sub2ind( size( Awv_if ), Aifridx, Acidx );

Awv_if( linifdx ) = Fif;

bwv_if = zeros( numWv, 1 );

end