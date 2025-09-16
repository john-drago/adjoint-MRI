function dGzdt = dGzdt_fn( t, s )
t = t( : ).';

dGzdt = d2kzdt2_fn( t, s ) / s.gyro;
end