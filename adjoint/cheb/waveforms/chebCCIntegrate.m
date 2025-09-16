function I = chebCCIntegrate( f, dom )
% This function will perform Clenshaw-Curtis quadrature of a function that
% is sampled at appropriate Chebyshev points of the second kind.

if nargin < 2
    dom = [ -1, 1 ];
end

n = size( f, 1 );
w = chebCCQuadWeights( n );

Iunsc = transpose( w ) * f;

I = ( diff( dom ) / 2 ) * Iunsc;

end