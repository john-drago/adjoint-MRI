function D = getChebDerivativeMatrix( numC, dom )
% This function will return the derivative matrix for the Chebyshev
% coefficients. 
%
% For a function:
% f(x) = c_0 T_0(x) + c_1 T_1(x) + c_2 T_2( x ) + ... + c_N T_N(x)
%
% This script will calculate the new coefficients for a Chebyshev
% polynomial based on the original function coefficients.
% df(x)/dx = c'_0 T_0 + c'_1 T_1( x ) + c'_2 T_2( x ) + ... + c'_{N-1} T_{N-1}(x)
% 
% c' = D c

D = zeros( (numC-1), numC );

for cc = (numC-1):-2:1
    D( cc  , (numC  ):-2:cc ) = 2 * ( (numC-1):-2:(cc  ) );
    if cc > 1
        D( cc-1, (numC-1):-2:cc ) = 2 * ( (numC-2):-2:(cc-1) );
    end
end

% Scale result for c'_0
D( 1, : ) = D( 1, : )/2;

if nargin > 1 % scale with domain
    dd = diff( dom );
    if dd <= 0
        error( "Needs to be increasing order in domain." )
    end
    % x( y ) = 2/(b-a) * ( y-a ) - 1
    % dx/dy = 2/(b-a);
    scDeriv = ( 2 / dd );
    D = scDeriv * D;
end

end