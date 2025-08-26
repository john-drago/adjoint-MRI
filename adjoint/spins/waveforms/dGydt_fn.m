function dGydt = dGydt_fn( t, s )
t = t( : ).';

dGydt = d2kydt2_fn( t, s ) / s.gyro;
end