function [ f, x ] = chebMinMax( c, dom, maxExtrema )
% Calculate global minimum and maximum of a Chebyshev representation of a
% polynomial represented by the coefficients in c. Function will calculate
% the derivative of the coefficients in c and then find the roots of this
% new Chebyshev polynomial. The function will then be evaluated at the
% roots of the derivative Chebyshev polynomial to determine the min and the
% max.
% 
% f(1) is minimum value and f(2) is maximum value
% the polynomial evaluated at x(1) gives f(1) and the polynomial evaluated
% at x(2) gives f(2).
%
% Also have the option to find numExtrema mins and maxes
% 
% Inspiration from CHEBTECH/MINANDMAX

if nargin < 2
    dom = [ -1 1 ];
end
if nargin < 3
    maxExtrema = 1;
end

[ f, x ] = findMinMax( c, dom, maxExtrema );

end

%% Helper Function
% ----------------------------------------------------------------------- %
function [ f, x ] = findMinMax( c, dom, maxExtrema )

nw = size( c, 2 );
if maxExtrema == 1
    f = zeros( 2, nw );
    x = zeros( 2, nw );
    oneExtrema = true;
elseif maxExtrema > 1
    f = repmat( {zeros( maxExtrema, 1 )}, [ 2, nw ] );
    x = repmat( {zeros( maxExtrema, 1 )}, [ 2, nw ] );
    oneExtrema = false;
else
    error( "Don't know how to process this number of extrema." )
end

cp = derChebCoeff( c, dom );

r = chebRoots( cp, dom );

for ww = 1:nw
    if ~isempty( r )
        if iscell( r )
            ri = [ dom(1); r{ ww }; dom(2) ];
        else
            ri = [ dom(1); r; dom(2) ];
        end
    else
        ri = [ dom(1); dom(2) ];
    end
    vi = evalChebClenshaw( ri, c( :, ww ), dom );
    
    % Assign minimums
    if oneExtrema
        [ f( 1, ww ), minidx ] = min( vi );
        x( 1, ww ) = ri( minidx );
    else
        [ minkVals, minIdxs ] = mink( vi, maxExtrema );
        num_mink = length( minIdxs );

        f{ 1, ww }( 1:num_mink ) = minkVals;
        x{ 1, ww }( 1:num_mink ) = ri( minIdxs );
    end

    % Assign maximums
    if oneExtrema
        [ f( 2, ww ), maxidx ] = max( vi );
        x( 2, ww ) = ri( maxidx );
    else
        [ maxkVals, maxIdxs ] = maxk( vi, maxExtrema );
        num_maxk = length( maxIdxs );

        f{ 2, ww }( 1:num_maxk ) = maxkVals;
        x{ 2, ww }( 1:num_maxk ) = ri( maxIdxs );
    end
end

end
% ----------------------------------------------------------------------- %

% % ----------------------------------------------------------------------- %
% function [ fi, xi ] = findMinMax( c, dom )
% 
% nw = size( c, 2 );
% 
% cpi = derChebCoeff( c, dom );
% 
% f = zeros( 2, nw );
% x = zeros( 2, nw );
% 
% r = chebRoots( cpi, dom );
% r = [ dom(1); r; dom(2) ];
% v = evalChebClenshaw( r, c, dom );
% 
% % Assign minimums
% [ fi( 1 ), minidx ] = min( v );
% xi( 1 ) = r( minidx );
% 
% % Assign maximums
% [ fi( 2 ), maxidx ] = max( v );
% xi( 2 ) = r( maxidx );
% 
% end
% % ----------------------------------------------------------------------- %