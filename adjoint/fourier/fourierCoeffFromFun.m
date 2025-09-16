function [ c, fvec ] = fourierCoeffFromFun( wvformsampfun, tdom, ncoeff )

pulseLength = diff( sort( tdom ) );
dt = pulseLength / ncoeff;
tsample = transpose( dt/2 : dt : pulseLength );

wvsample = wvformsampfun( tsample );

c = fftshift( fft( wvsample, ncoeff, 1 ), 1 );
c = c / ncoeff; % for scaling purposes although technically not correct

shift = true;
fvec = transpose( fvecDFT( 1/dt, ncoeff, shift ) );

end