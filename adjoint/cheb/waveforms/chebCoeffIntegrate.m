function I = chebCoeffIntegrate( c, dom )
% This function will integrate a chebyshev function represented by
% coefficients. Please see CHEBTECH/SUM of chefun for more details.

if nargin < 2
    dom = [ -1, 1 ];
end

n = size( c, 1 ); % number of chebyshev coefficients

% If only one input point
if ( n == 0 )
    error( "Don't know how to parse ." )
elseif ( n == 1 )
    Iunsc = 2 * c;
else
    % See Theorem 19.2 of ATAP for explanation
    c( (2:2:end), :) = 0;
    Iunsc = [ 2, 0, 2./( 1 - (2:(n-1)).^2) ] * c;
end

I = ( diff( dom ) / 2 ) * Iunsc;
end