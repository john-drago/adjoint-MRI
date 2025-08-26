function cpp = der2ChebCoeff( c, dom )
% Function will calculate the second derivative of a function represented
% by chebyshev coefficients. Can evaluate multiple waveforms that are
% represented by a column in c.
%
% For a function:
% f(x) = c_0 T_0(x) + c_1 T_1(x) + c_2 T_2( x ) + ... + c_N T_N(x)
%
% This script will calculate the new coefficients for a Chebyshev
% polynomial based on the original function coefficients.
% d^2f(x)/dx^2 = c''_0 T_0 + c''_1 T_1( x ) + c''_2 T_2( x ) + ... + c''_{N-2} T_{N-2}(x)
%
% This is performed by using the following relation:
%
% c''_{i-1} = c''_{i+1} + 2 i c'_i, for i = N-1,...,i
% with initial conditions: c''_{i+1} = 0 and c''_{i+2} = 0
%
% Derived in Section 2.4.4 and 2.4.5 of Chebyshev Polynomials by J.C. Mason
%
% This can be expanded to:
%
% For i odd >= 1:
% c'_i = \sum^N_{ i even} \sum^N_{ j > i odd } 2*( i ) * 2 *( j ) c_j
%
% For i even > 1:
%
% c'_i = \sum^N_{ i odd} \sum^N_{ j > i even } 2*( i ) * 2 *( j ) c_j
%
% For i = 0:
% c''_0 = c''_2/2 + c'_1
%
% Function is inspired by chebtech2/diff function of chebfun.

[ ne_init, nw ] = size( c );
nc_init = ne_init - 1;
nc = nc_init - 2;
ne = nc + 1;

cpp = zeros( ne, nw );
w = repmat( 2 * ( (1:(nc_init)).'), [ 1, nw ] );
wc = zeros( ne, nw );

wc( (end  ):-2:1, : ) = w( (end-1):-2:1, : ) .* cumsum( w( (end  ):-2:2, : ) .* c( (end  ):-2:3, : ) ); % weighting for c_{ N   }, c_{ N-2 }, c_{ N-4 }
wc( (end-1):-2:1, : ) = w( (end-2):-2:2, : ) .* cumsum( w( (end-1):-2:2, : ) .* c( (end-1):-2:3, : ) ); % weighting for c_{ N-1 }, c_{ N-3 }, c_{ N-5 }

cpp( (end  ):-2:1, : ) = cumsum( wc( (end  ):-2:1, : ), 1 );
cpp( (end-1):-2:1, : ) = cumsum( wc( (end-1):-2:1, : ), 1 );

% rescale c''_0
cpp( 1, : ) = 0.5 * cpp( 1, : );

if nargin > 1 % scale with domain
    dd = diff( dom );
    if dd <= 0
        error( "Needs to be increasing order in domain." )
    end
    % x( y ) = 2/(b-a) * ( y-a ) - 1
    % dx/dy = 2/(b-a);
    scDeriv = ( 2 / dd )^2;
    cpp = scDeriv * cpp;
end

end