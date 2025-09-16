function Gx = Gx_fn( t, s )
Gx = dkxdt_fn( t, s ) / s.gyro;
end