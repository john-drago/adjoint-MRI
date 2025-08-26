function f = evalChebClenshaw( x, c, dom )
% Function will evaluate chebyshev coefficients at points in x. Can
% evaluate multiple waveforms that have coefficients in columns of c
%
% Below code is heavily inspired by Nick Hale's implementation in chebfun
% function chebtech/CLENSHAW written. Specifically, we are doing two iterations of
% the algortihm at once.
%
% For a function:
%
% f(x) = c_0 T_0(x) + c_1 T_1(x) + c_2 T_2( x ) + ... + c_N T_N(x)
%
% This can be evaluated using Clenshaw's algorithm.
%
% alpha(x) = 2x, beta = -1
%
% b_k(x) = c_k + 2x b_{k+1}(x) - b_{k+2}(x)
%
% f(x) = c_0 + x b_1(x) - b_2(x)

x = x( : ); % ensure column vector

if nargin > 2 % scale with domain
    dd = diff( dom );
    if dd <= 0
        error( "Needs to be increasing order in domain." )
    end
    x = chebFwdMap( x, dom );
end

nt = size( x, 1 );
nc = size( c, 1 );
nw = size( c, 2 );

x = repmat( x, [ 1, nw ]);
bk1 = zeros( nt, nw );
bk2 = zeros( nt, nw );
twox = 2*x;
e = ones(size(x, 1), 1);

for kk = ( (nc):-2:3 )
    bk2 = e * c(kk,:) + twox.*bk1 - bk2;
    bk1 = e * c(kk-1,:) + twox.*bk2 - bk1;
end
if ( mod((nc-1), 2) )
    tmp = bk1;
    bk1 = e*c(2,:) + twox.*bk1 - bk2;
    bk2 = tmp;
end
f = e * c(1,:) + x.*bk1 - bk2;

end