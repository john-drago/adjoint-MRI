function DFB = getFourierBasisDerivativeMatrix( orderFourier, tvec, omega0 )

orderFourier = double( orderFourier );

k = 0 : ( orderFourier );

sinkOmega = sin( tvec( : ) * ( k * omega0 ) );
coskOmega = cos( tvec( : ) * ( k( 2:end ) * omega0 ) );

dk = [ -( 0 : orderFourier ), ( 1 : orderFourier ) ]  ;
DFB = ( dk * omega0 ) .* [ sinkOmega( :, 1:( orderFourier + 1 ) ), coskOmega( :, 1:( orderFourier ) ) ];

end