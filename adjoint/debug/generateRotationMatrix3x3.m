function R = generateRotationMatrix3x3( B, opt )
dt = opt.dt;
gyro = opt.gyro;

Bmag = sqrt( B( 1 ).^2 + B( 2 ).^2 + B( 3 ).^2 );

u = B ./ Bmag;
u= u(:);

phi = -gyro * dt * Bmag;

R = cos( phi ) * eye( 3 ) + ( 1 - cos( phi ) ) * (u * u.') + sin( phi ) *...
    generateCrossProductMatrix( u );

end