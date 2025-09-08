function [ ax, ang ] = determineAxisAngleFromRotMat( R )
ax = null( R - eye( 3 ) );
V = generateCrossProductMatrix( ax );
cosang = ( trace( R ) - 1 ) / 2;
sinang = - ( trace( V * R ) / 2 );
ang = atan2( sinang, cosang );
end