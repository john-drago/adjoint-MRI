function D = getFourierRFBasisDerivativeMatrix( orderFourier, FB, omega0 )

orderFourier = double( orderFourier );

k = ( -orderFourier : 1 : orderFourier );
D = 1j * ( ( k * omega0 ) ) .* FB;

end