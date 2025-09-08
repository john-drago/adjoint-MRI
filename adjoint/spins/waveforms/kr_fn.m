function kr = kr_fn( t, s )
t = t( : ).';

kr = s.kmax ./ denom_fn( t, s );
end