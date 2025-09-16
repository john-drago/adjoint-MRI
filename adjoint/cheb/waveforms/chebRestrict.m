function cr = chebRestrict( c, s, dom, absflag )
% Function will restrict a function expressed in chebyshev coefficients to
% a subinterval specified in s of the domain specified in dom.
%
% CHEBTECH/RESTRICT

% check domain
domtol = 1e-12;
s = s( : );
if any( s < ( dom(1) - domtol ) ) || any( s > ( dom(end) + domtol ) )
    error( "Restriction points outside the bounds of the domain." )
end

% initialize cell array
[ nc, nw ] = size( c );
ns = length( s ) - 1;
cr = cell( ns, nw );

% convert points of s to [ -1 1 ] domain
st = chebFwdMap( s, dom );

% Sample nc points in [-1 1 ] domain
x = chebPts2( nc );

% now, determine locations of nc points in each of the domains specified by
% s 
y = 0.5 * [ 1 - x, 1 + x ] * transpose( [ st(1:(end-1), 1 ), st(2:(end), 1 ) ] );
yvec = reshape( y, [ nc*ns, 1 ] );

for ww = 1:nw
    % evaluate polynomial at these points in each of the subdomains
    valvec = evalChebClenshaw( yvec, c( :, ww ) );
    val = reshape( valvec, [ nc, ns ] );
    if absflag
        val = abs( val );
    end

    % convert each of these subdomains to their own chebyshev coefficients
    cw = chebValToCoeff( val );
    cr( :, ww ) = transpose( mat2cell( cw , nc, ones( 1, ns ) ) );
end

end