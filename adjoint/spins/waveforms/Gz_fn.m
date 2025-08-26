function Gz = Gz_fn( t, s )
Gz = dkzdt_fn( t, s ) / s.gyro;
end