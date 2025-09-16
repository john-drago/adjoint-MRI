function DFB = getFourierBasis2DerivativeMatrix( orderFourier, tvec, omega0 )

orderFourier = double( orderFourier );

k = 0 : ( orderFourier );

coskOmega = cos( tvec( : ) * ( k * omega0 ) );
sinkOmega = sin( tvec( : ) * ( k( 2:end ) * omega0 ) );

d2k = -[ ( ( 0 : orderFourier ).^2 ), ( 1 : orderFourier ).^2 ]  ;
DFB = ( d2k * ( omega0^2 ) ) .* [ coskOmega( :, 1:( orderFourier + 1 ) ), sinkOmega( :, 1:( orderFourier ) ) ];

end