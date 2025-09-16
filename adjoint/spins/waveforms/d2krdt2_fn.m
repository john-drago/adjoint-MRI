function d2krdt2 = d2krdt2_fn( t, s )
t = t( : ).';

denom = denom_fn( t, s );
indenom = indenom_fn( t, s );

d2krdt2 = ( ( s.a.^2 * s.kmax * indenom ) .* ( indenom - 1 ) ) ./ ( s.T^2 * denom.^3 );
end