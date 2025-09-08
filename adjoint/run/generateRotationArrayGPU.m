function Rarray = generateRotationArrayGPU( Rarray, wv )
% This function will generate the rotation array 
% Rarray dims: Nvox x Nt x 9 

Bx = wv.b1preal * wv.breal - wv.b1pimag * wv.bimag;
By = wv.b1pimag * wv.breal + wv.b1preal * wv.bimag;

if wv.numZCoils > 0
    if ~wv.constantRotatingFrame
        Bz = wv.bzsens * wv.Shim + wv.pos * wv.Grad + wv.db0 + (1/wv.gyro) * wv.dwxyvec;
    else
        Bz = wv.bzsens * wv.Shim + wv.pos * wv.Grad + wv.db0;
    end
else
    if ~wv.constantRotatingFrame
        Bz = wv.pos * wv.Grad + wv.db0 + (1/wv.gyro) * wv.dwxyvec;
    else
        Bz = wv.pos * wv.Grad + wv.db0;
    end
end

magB = sqrt( Bx.^2 + By.^2 + Bz.^2 );
phi = - magB .* ( wv.gyro * wv.dtvec );

ux = Bx ./ magB;
uy = By ./ magB;
uz = Bz ./ magB;

ux( magB < 1e-12 ) = 0;
uy( magB < 1e-12 ) = 0;
uz( magB < 1e-12 ) = 0;

uxuy = ux .* uy;
uxuz = ux .* uz;
uyuz = uy .* uz;

uxsq = ux.^2;
uysq = uy.^2;
uzsq = uz.^2;

cosphi = cos( phi );
onemcosphi = ( 1 - cosphi );

sinphi = sin( phi );
uxsinphi = sinphi .* ux;
uysinphi = sinphi .* uy;
uzsinphi = sinphi .* uz;
uxuyonemcosphi = uxuy .* onemcosphi;
uxuzonemcosphi = uxuz .* onemcosphi;
uyuzonemcosphi = uyuz .* onemcosphi;

Rarray( :, :, 1 ) = cosphi + onemcosphi .* uxsq;
Rarray( :, :, 2 ) = uxuyonemcosphi + uzsinphi;
Rarray( :, :, 3 ) = uxuzonemcosphi - uysinphi;
Rarray( :, :, 4 ) = uxuyonemcosphi - uzsinphi;
Rarray( :, :, 5 ) = cosphi + onemcosphi .* uysq;
Rarray( :, :, 6 ) = uyuzonemcosphi + uxsinphi;
Rarray( :, :, 7 ) = uxuzonemcosphi + uysinphi;
Rarray( :, :, 8 ) = uyuzonemcosphi - uxsinphi;
Rarray( :, :, 9 ) = cosphi + onemcosphi .* uzsq;


end