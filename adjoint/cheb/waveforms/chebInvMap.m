function y = chebInvMap( x, dom )
% Map from x in [ -1 1 ] to y in [ a b ]
dd = diff( dom );
if dd <= 0
    error( "Needs to be increasing order in domain." )
end
y = ( x .* dd + sum( dom ) ) ./ 2;

y = max( y, dom( 1 ) );
y = min( y, dom( 2 ) );

end