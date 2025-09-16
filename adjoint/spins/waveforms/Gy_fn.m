function Gy = Gy_fn( t, s )
Gy = dkydt_fn( t, s ) / s.gyro;
end