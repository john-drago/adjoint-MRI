function cullIdxLin = cullSumSquaresCorrelationSpokes(...
    p, resid, b1pnormalized, spokeIdx, Kx, Ky, lini_sort, opt, spokes )

numIdx = size( lini_sort, 1 );
corrList = zeros( numIdx, 1 );

for ll = 1:numIdx
    ki = [ Kx( lini_sort( ll ) ); Ky( lini_sort( ll ) ) ];
    lcols = makeSpokesSTAColumns( b1pnormalized, ki, spokeIdx, opt, spokes );
    innprod = ctranspose( resid ) * lcols;
    corrList( ll ) = sum( abs(innprod).^2 );
end

[ ~, corrListIdx ] = sort( corrList, 'descend' );

cullIdxLin = lini_sort( corrListIdx( 1:p ) );

end