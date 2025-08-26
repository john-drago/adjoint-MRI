function [ FB, fvec ] = getFourierBasisMatrix( orderFourier, tvec, omega0 )

orderFourier = double( orderFourier );

f0 = omega0 / ( 2*pi );

k = 0 : orderFourier;

coskOmega = cos( tvec( : ) * ( k * omega0 ) );
sinkOmega = sin( tvec( : ) * ( k( 2:end ) * omega0 ) );

FB = [ coskOmega( :, 1:( orderFourier + 1 ) ), sinkOmega( :, 1:( orderFourier ) ) ];
fvec = transpose(...
    [ ( k( 1:( orderFourier + 1 ) ) * f0 ), ( k( 2:( orderFourier + 1 ) ) * f0 ) ]...
    );

end