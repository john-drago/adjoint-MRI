function w = chebCCQuadWeights( n )
% get quadrature weights at the Chebyshev points of the second kind for
% Clenshaw-Curtis quadrature. Code inspiration taken directly from 
% CHEBTECH2/QUADWTS of chebfun. Will return as a column vector.
% 
% Please see Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
% quadrature rules", 2006 for further information and von Winckel's file
% exchange: 
% https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature

if ( n == 0 ) % Deal with no points
    error( "Don't know how to parse no points for Clenshaw-Curtis quadrature." )
elseif ( n == 1 ) % Deal with single point
    w = 2;
    
else % Deal with other cases

    % Exact integrals of T_k (even)
    c = 2./[1;  1-transpose( (2:2:(n-1)).^2 )];

    % Mirror for DCT via FFT
    c = [c; c(floor(n/2):-1:2)];

    % Interior weights
    w = ifft(c);

    % adjust for the boundary points
    w( [ 1, n ] ) = w( 1 ) / 2;
end

end