function ktheta = ktheta_fn( t, s )
t = t( : ).';

ktheta = s.u * t;
end