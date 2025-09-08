function dGxyzdt = dGxyzdt_fn( t, s )
t = t( : ).';

dGxyzdt = d2kxyzdt2_fn( t, s ) / s.gyro;
end