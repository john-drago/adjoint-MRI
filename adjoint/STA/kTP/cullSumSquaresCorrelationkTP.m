function cullIdxLin = cullSumSquaresCorrelationkTP( p, resid, b1pnormalized, kTTime, Kx, Ky, Kz, lini_sort, opt, pulse )

numIdx = size( lini_sort, 1 );
corrList = zeros( numIdx, 1 );

for ll = 1:numIdx
    ki = [ Kx( lini_sort( ll ) ); Ky( lini_sort( ll ) ); Kz( lini_sort( ll ) ) ];
    lcols = makekTPSTAColumns( b1pnormalized, ki, kTTime, opt, pulse );
    innprod = resid' * lcols;
    corrList( ll ) = sum( abs(innprod).^2 );
end

[ ~, corrListIdx ] = sort( corrList, 'descend' );

cullIdxLin = lini_sort( corrListIdx( 1:p ) );

end