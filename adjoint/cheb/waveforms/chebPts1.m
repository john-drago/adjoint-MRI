function x = chebPts1( N )
% Function that will return N chebpoints of the first kind.
%
% See CHEBTECH1/CHEBPTS from chebfun.
% 
% See Chapter 2 of ATAP

% Use sine for symmetry
if ( N == 0 ) % Special case (no points)
    error( "Don't know how to process no points" )
elseif ( N == 1 ) % single point
    x = 0;
else
    N = double(N);
    x = transpose( sin( pi * ( ( (-N+1) : 2 : (N-1) ) / ( 2 * N ) ) ) );
end
end