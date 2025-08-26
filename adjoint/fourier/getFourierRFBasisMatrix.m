function [ FB, fvec ] = getFourierRFBasisMatrix( orderFourier, tvec, omega0 )

f0 = omega0 / ( 2*pi );

orderFourier = double( orderFourier );

k = -orderFourier : 1 : orderFourier;
expjkOmega = exp( 1j * tvec( : ) * ( k * omega0 ) );

FB = expjkOmega;
fvec = transpose( ( k * f0 ) );

end