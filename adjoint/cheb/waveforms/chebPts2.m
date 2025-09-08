function x = chebPts2( N )
% Function that will return N chebpoints of the second kind.
%
% See CHEBTECH2/CHEBPTS from chebfun.
% 
% See Chapter 2 of ATAP for an exercise to explain this implementation

% Use sine for symmetry

if ( N == 0 ) % Special case (no points)
    error( "Don't know how to process no points" )
elseif ( N == 1 ) % single point
    x = 0;
else
    M = double(N)-1;
    x = transpose( sin( pi * ( (-M) : 2 : M ) / ( 2 * M ) ) );
end


end