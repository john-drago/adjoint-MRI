function R = generateRotationMatrixAxisAngle( ax, ang )

ax = ax( : );
cosang = cos( ang );
sinang = sin( ang );

R =...
    cosang * eye( 3 ) +...
    sinang * generateCrossProductMatrix( ax ) +...
    ( ( 1 - cosang ) ) *  ( ax * ax' );

end