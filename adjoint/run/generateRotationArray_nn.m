function [ R_nn, B ] = generateRotationArray_nn( B, R_nn, wv, nn )

B( :, 1 ) = wv.b1preal * wv.breal( :, nn ) - wv.b1pimag * wv.bimag( :, nn );
B( :, 2 ) = wv.b1pimag * wv.breal( :, nn ) + wv.b1preal * wv.bimag( :, nn );

if ( wv.numZCoils > 0 )
    if ~wv.constantRotatingFrame && ( wv.dwxyvec( nn ) ~= 0 )
        B( :, 3 ) = wv.bzsens * wv.Shim( :, nn ) + wv.pos * wv.Grad( :, nn ) + ( wv.db0 + ( 1/wv.gyro ) * wv.dwxyvec( nn ) );
    else
        B( :, 3 ) = wv.bzsens * wv.Shim( :, nn ) + wv.pos * wv.Grad( :, nn ) + wv.db0;
    end
else
    if ~wv.constantRotatingFrame && ( wv.dwxyvec( nn ) ~= 0 )
        B( :, 3 ) = wv.pos * wv.Grad( :, nn ) + ( wv.db0 + ( 1/wv.gyro ) * wv.dwxyvec( nn ) );
    else
        B( :, 3 ) = wv.pos * wv.Grad( :, nn ) + wv.db0;
    end
end

magB_nn = sqrt( B( :, 1 ).^2 + B( :, 2 ).^2 + B( :, 3 ).^2 );

phi = - ( wv.gyro * wv.dtvec( nn ) ) * magB_nn; 
u = B ./ magB_nn ;
u( magB_nn < 1e-12 ) = 0;

% Precompute values
uxuy = u( :, 1 ) .* u( :, 2 );
uxuz = u( :, 1 ) .* u( :, 3 );
uyuz = u( :, 2 ) .* u( :, 3 );
usq = u.^2;
cosphi = cos( phi );
onemcosphi = ( 1 - cosphi );

sinphi = sin( phi );
usinphi = sinphi .* u;
uxuyonemcosphi = uxuy .* onemcosphi;
uxuzonemcosphi = uxuz .* onemcosphi;
uyuzonemcosphi = uyuz .* onemcosphi;

R_nn( :, 1 ) = cosphi + onemcosphi .* usq( :, 1 );
R_nn( :, 2 ) = uxuyonemcosphi + usinphi( :, 3 );
R_nn( :, 3 ) = uxuzonemcosphi - usinphi( :, 2 );
R_nn( :, 4 ) = uxuyonemcosphi - usinphi( :, 3 );
R_nn( :, 5 ) = cosphi + onemcosphi .* usq( :, 2 );
R_nn( :, 6 ) = uyuzonemcosphi + usinphi( :, 1 );
R_nn( :, 7 ) = uxuzonemcosphi + usinphi( :, 2 );
R_nn( :, 8 ) = uyuzonemcosphi - usinphi( :, 1 );
R_nn( :, 9 ) = cosphi + onemcosphi .* usq( :, 3 );

end