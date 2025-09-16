function x = chebFwdMap( y, dom )
% Map from y in [ a b ] to x in [ -1 1 ]
dd = diff( dom );
if dd <= 0
    error( "Needs to be increasing order in domain." )
end
x = ( 2 .* y - sum( dom ) ) ./ dd;
end