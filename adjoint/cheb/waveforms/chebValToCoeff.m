function c = chebValToCoeff( val )
% This function will take the values of a function samples at chebyshev
% points of the second kind and return the coefficients in an expansion of
% chebyshev polynomials of the first kind.
%
% This function was inspired in part by Chebfun's chebtech2/coeffs2vals
% function
%
% For discussion of how the FFT can be incorporated to determine how
% coefficients can be determined from samples of cheb polynomials of the
% first kind, please see "Interpolation by Fast Fourier and Chebyshev
% Transforms", International Journal for Numerical Methods in Engineering,
% Monro, 1979.
%
% Equations 27 and 28 Monro 
% For F_k = real( (1/2) \sum_{n=0}^{2N-1} C_n e^{ j 2 pi k n / (2N) } ), k = 0,1,..., 2N-1
%     C_n = real( 2/(2N) \sum_{k=0}^{2N-1} F_k e^{ -j 2 pi k n/(2N) } ), n = 0,1,..., 2N-1
% 

[ valn, valw ] = size( val );

% Flip the function values, because we follow the Chebfun standard, where
% we sample Cheb points of the second kind in order:
%
% xk = cos( k*pi/(N) ), k = N,N-1,...,0
valexp = zeros( (2*valn-2), valw );
valexp( 1:valn, : ) = val( valn:-1:1, : );
valexp( (valn+1):(2*valn-2), : ) = val( 2:1:(valn-1), : );

if ( isreal(valexp) ) % real values
    % Real-valued case:
    c = ifft( valexp );
    c = real( c );
elseif ( isreal(1j * valexp) ) % imaginary values
    c = ifft( imag( valexp ) );
    c = 1j * real( c );
else % could be complex
    c = ifft( valexp );
end

c = c( 1:valn, : ); % get first N values
c( 2:(valn-1), : ) = 2 * c( 2:(valn-1), : ); % scale. See Equation 31 Monro

end