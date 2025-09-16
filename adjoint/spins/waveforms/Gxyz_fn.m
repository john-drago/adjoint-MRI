function Gxyz = Gxyz_fn( t, s )
t = t( : ).';
Gxyz = dkxyzdt_fn( t, s ) / s.gyro;
end