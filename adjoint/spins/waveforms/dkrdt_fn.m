function dkrdt = dkrdt_fn( t, s )
t = t( : ).';
dkrdt = ( ( -s.kmax * s.a ) * indenom_fn( t, s ) ) ./ ( s.T * denom_fn( t, s ).^2 );
end