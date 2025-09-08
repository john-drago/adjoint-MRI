function cp = derChebCoeff( c, dom )
% Function will calculate the derivative of a function represented by chebyshev
% coefficients. Can evaluate multiple waveforms that are represented by a
% column in c.
%
% For a function:
% f(x) = c_0 T_0(x) + c_1 T_1(x) + c_2 T_2( x ) + ... + c_N T_N(x)
%
% This script will calculate the new coefficients for a Chebyshev
% polynomial based on the original function coefficients.
% df(x)/dx = c'_0 T_0 + c'_1 T_1( x ) + c'_2 T_2( x ) + ... + c'_{N-1} T_{N-1}(x)
%
% This is performed by using the following relation:
%
% c'_{i-1} = c'_{i+1} + 2 i c_i, for i = N-1,...,i
% with initial conditions: c'_{i+1} = 0 and c'_{i+2} = 0
%
% Derived in Section 2.4.4 and 2.4.5 of Chebyshev Polynomials by J.C. Mason
%
% This can be expanded to:
%
% For i odd >= 1:
% c'_i = \sum^N_{ i even} 2*( i ) c_i
%
% For i even > 1:
% c'_i = \sum^N_{ i odd} 2*( i ) c_i
%
% 
% For i = 0:
% c'_0 = c'_2/2 + c_1
%
% Function is inspired by chebtech2/diff function of chebfun.


[ nc, nw ] = size( c );

if nc > 1
    ncp = nc - 1; % num coefficients for cp
    cp = zeros( ncp, nw );
    w = repmat( ( 2 * (1:(ncp)).'), [ 1, nw ] );
    wc = w .* c( 2:end, : );
    cp( (ncp):-2:1, : ) = cumsum( wc( (ncp):-2:1, : ), 1 );
    cp( (ncp-1):-2:1, : ) = cumsum( wc( (ncp-1):-2:1, : ), 1 );
    cp( 1, : ) = 0.5 * cp( 1, : );

    if nargin > 1 % scale with domain
        dd = diff( dom );
        if dd <= 0
            error( "Needs to be increasing order in domain." )
        end
        % x( y ) = 2/(b-a) * ( y-a ) - 1
        % dx/dy = 2/(b-a);
        scDeriv = ( 2 / dd );
        cp = scDeriv * cp;
    end
else
    cp = zeros( 1, nw );
end

end