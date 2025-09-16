function val = chebCoeffToVal( c )
% This function will take the chebyshev coefficient values for chebyshev
% expansion in chebyshev polynomials of the first kind and will return
% evaluations at the cheb points on the second kind.
%
% This function was inspired in part by Chebfun's chebtech2/coeffs2vals
% function
%
% For discussion of how the FFT can be incorporated to evaluate the
% chebyshev polynomials given their coefficients, please see 
% "Interpolation by Fast Fourier and Chebyshev Transforms", International
% Journal for Numerical Methods in Engineering, Monro, 1979.
%
% Equations 27 and 28 Monro
% For F_k = real( (1/2) \sum_{n=0}^{2N-1} C_n e^{ j 2 pi k n / (2N) } ), k = 0,1,..., 2N-1
%     C_n = real( 2/(2N) \sum_{k=0}^{2N-1} F_k e^{ -j 2 pi k n/(2N) } ), n = 0,1,..., 2N-1
% 
% Assume domain is from -1 to 1 or otherwise defined by domain = [ a, b ];

[ cn, cw ] = size( c );

c( 2:(cn-1), : ) = c( 2:(cn-1), : )/2; % scaling. See Equation 30 in Monro

cexp = zeros( 2*cn-2, cw ); 
cexp( 1:cn, : ) = c;
cexp( (cn+1):(2*cn-2), : ) = c( (cn-1):-1:2, : );

if ( isreal(c) ) % real values
    val = real(fft( cexp ));
elseif ( isreal( 1j * cexp ) ) % imaginary values
    val = 1j*real(fft(imag( cexp )));
else % could be complex
    val = fft( cexp );
end

% Flip the function values, because we follow the Chebfun standard, where
% we sample Cheb points of the second kind in order:
%
% xk = cos( k*pi/(N) ), k = N,N-1,...,0
val = val( (cn):-1:1, : );

end