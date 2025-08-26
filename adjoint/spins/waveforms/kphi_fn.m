function kphi = kphi_fn( t, s )
t = t( : ).';

kphi = s.v * t;
end