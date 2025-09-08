function dGxdt = dGxdt_fn( t, s )
t = t( : ).';

dGxdt = d2kxdt2_fn( t, s ) / s.gyro;
end