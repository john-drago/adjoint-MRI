function initSt = calcpartialfpartialBRotationMatrixCPU( initSt, wv, Mnn, Mnnp1, nn )
% This function will calculate the partialfpartialB assuming the use of a
% rotation matrix for forward model integration.


%% use cell array approach
Mttx = Mnn( :, 1 );
Mtty = Mnn( :, 2 );
Mttz = Mnn( :, 3 );

Mttp1x = Mnnp1( :, 1 );
Mttp1y = Mnnp1( :, 2 );
Mttp1z = Mnnp1( :, 3 );

initSt.B{ 1 } = wv.b1preal * wv.breal( :, nn ) - wv.b1pimag * wv.bimag( :, nn );
initSt.B{ 2 } = wv.b1pimag * wv.breal( :, nn ) + wv.b1preal * wv.bimag( :, nn );
initSt.B{ 3 } = wv.bzsens * wv.Shim( :, nn ) + wv.pos * wv.Grad( :, nn ) + wv.db0;

magBtt = sqrt( initSt.B{ 1 }.^2 + initSt.B{ 2 }.^2 + initSt.B{ 3 }.^2 );

dt = wv.dtvec( nn );

phi = - ( wv.gyro * dt ) * magBtt; 
ux = initSt.B{ 1 } ./ magBtt ;
uy = initSt.B{ 2 } ./ magBtt ;
uz = initSt.B{ 3 } ./ magBtt ;

% Precompute values
uxuy = ux .* uy;
uxuz = ux .* uz;
uyuz = uy .* uz;

uxsq = ux.^2;
uysq = uy.^2;
uzsq = uz.^2;

cosphi = cos( phi );
onemcosphi = ( 1 - cosphi );

sinphi = sin( phi );

udotM = ux .* Mttx + uy .* Mtty + uz .* Mttz;

% start assembling gradients
initSt.partialfpartialu{ 1 } = onemcosphi .* ( udotM + ux .* Mttx );
initSt.partialfpartialu{ 4 } = onemcosphi .* ( ux .* Mtty ) + sinphi .* Mttz;
initSt.partialfpartialu{ 7 } = onemcosphi .* ( ux .* Mttz ) - sinphi .* Mtty;
initSt.partialfpartialu{ 2 } = onemcosphi .* ( uy .* Mttx ) - sinphi .* Mttz;
initSt.partialfpartialu{ 5 } = onemcosphi .* ( udotM + uy .* Mtty );
initSt.partialfpartialu{ 8 } = onemcosphi .* ( uy .* Mttz ) + sinphi .* Mttx;
initSt.partialfpartialu{ 3 } = onemcosphi .* ( uz .* Mttx ) + sinphi .* Mtty;
initSt.partialfpartialu{ 6 } = onemcosphi .* ( uz .* Mtty ) - sinphi .* Mttx;
initSt.partialfpartialu{ 9 } = onemcosphi .* ( udotM + uz .* Mttz );

initSt.partialfpartialphi{ 1 } = - uz .* Mttp1y + uy .* Mttp1z;
initSt.partialfpartialphi{ 2 } = + uz .* Mttp1x - ux .* Mttp1z;
initSt.partialfpartialphi{ 3 } = - uy .* Mttp1x + ux .* Mttp1y;

initSt.partialupartialB{ 1 } = ( 1 - uxsq ) ./  magBtt;
initSt.partialupartialB{ 4 } = ( - uxuy ) ./  magBtt;
initSt.partialupartialB{ 7 } = ( - uxuz ) ./  magBtt;
initSt.partialupartialB{ 5 } = ( 1 - uysq ) ./  magBtt;
initSt.partialupartialB{ 8 } = ( - uyuz ) ./  magBtt;
initSt.partialupartialB{ 9 } = ( 1 - uzsq ) ./  magBtt;
initSt.partialupartialB{ 2 } = initSt.partialupartialB{ 4 };
initSt.partialupartialB{ 3 } = initSt.partialupartialB{ 7 };
initSt.partialupartialB{ 6 } = initSt.partialupartialB{ 8 };

initSt.partialphipartialB{ 1 } = ( -wv.gyro * dt ) * ux;
initSt.partialphipartialB{ 2 } = ( -wv.gyro * dt ) * uy;
initSt.partialphipartialB{ 3 } = ( -wv.gyro * dt ) * uz;

initSt.partialfpartialB( :, 1 ) = ...
    initSt.partialfpartialu{ 1 } .* initSt.partialupartialB{ 1 } +...
    initSt.partialfpartialu{ 4 } .* initSt.partialupartialB{ 2 } +...
    initSt.partialfpartialu{ 7 } .* initSt.partialupartialB{ 3 } +...
    initSt.partialfpartialphi{ 1 } .* initSt.partialphipartialB{ 1 };
initSt.partialfpartialB( :, 4 ) = ...
    initSt.partialfpartialu{ 1 } .* initSt.partialupartialB{ 4 } +...
    initSt.partialfpartialu{ 4 } .* initSt.partialupartialB{ 5 } +...
    initSt.partialfpartialu{ 7 } .* initSt.partialupartialB{ 6 } +...
    initSt.partialfpartialphi{ 1 } .* initSt.partialphipartialB{ 2 };
initSt.partialfpartialB( :, 7 ) = ...
    initSt.partialfpartialu{ 1 } .* initSt.partialupartialB{ 7 } +...
    initSt.partialfpartialu{ 4 } .* initSt.partialupartialB{ 8 } +...
    initSt.partialfpartialu{ 7 } .* initSt.partialupartialB{ 9 } +...
    initSt.partialfpartialphi{ 1 } .* initSt.partialphipartialB{ 3 };
initSt.partialfpartialB( :, 2 ) = ...
    initSt.partialfpartialu{ 2 } .* initSt.partialupartialB{ 1 } +...
    initSt.partialfpartialu{ 5 } .* initSt.partialupartialB{ 2 } +...
    initSt.partialfpartialu{ 8 } .* initSt.partialupartialB{ 3 } +...
    initSt.partialfpartialphi{ 2 } .* initSt.partialphipartialB{ 1 };
initSt.partialfpartialB( :, 5 ) = ...
    initSt.partialfpartialu{ 2 } .* initSt.partialupartialB{ 4 } +...
    initSt.partialfpartialu{ 5 } .* initSt.partialupartialB{ 5 } +...
    initSt.partialfpartialu{ 8 } .* initSt.partialupartialB{ 6 } +...
    initSt.partialfpartialphi{ 2 } .* initSt.partialphipartialB{ 2 };
initSt.partialfpartialB( :, 8 ) = ...
    initSt.partialfpartialu{ 2 } .* initSt.partialupartialB{ 7 } +...
    initSt.partialfpartialu{ 5 } .* initSt.partialupartialB{ 8 } +...
    initSt.partialfpartialu{ 8 } .* initSt.partialupartialB{ 9 } +...
    initSt.partialfpartialphi{ 2 } .* initSt.partialphipartialB{ 3 };
initSt.partialfpartialB( :, 3 ) = ...
    initSt.partialfpartialu{ 3 } .* initSt.partialupartialB{ 1 } +...
    initSt.partialfpartialu{ 6 } .* initSt.partialupartialB{ 2 } +...
    initSt.partialfpartialu{ 9 } .* initSt.partialupartialB{ 3 } +...
    initSt.partialfpartialphi{ 3 } .* initSt.partialphipartialB{ 1 };
initSt.partialfpartialB( :, 6 ) = ...
    initSt.partialfpartialu{ 3 } .* initSt.partialupartialB{ 4 } +...
    initSt.partialfpartialu{ 6 } .* initSt.partialupartialB{ 5 } +...
    initSt.partialfpartialu{ 9 } .* initSt.partialupartialB{ 6 } +...
    initSt.partialfpartialphi{ 3 } .* initSt.partialphipartialB{ 2 };
initSt.partialfpartialB( :, 9 ) = ...
    initSt.partialfpartialu{ 3 } .* initSt.partialupartialB{ 7 } +...
    initSt.partialfpartialu{ 6 } .* initSt.partialupartialB{ 8 } +...
    initSt.partialfpartialu{ 9 } .* initSt.partialupartialB{ 9 } +...
    initSt.partialfpartialphi{ 3 } .* initSt.partialphipartialB{ 3 };


initSt.partialfpartialBx = initSt.partialfpartialB( :, 1:3 );
initSt.partialfpartialBy = initSt.partialfpartialB( :, 4:6 );
initSt.partialfpartialBz = initSt.partialfpartialB( :, 7:9 );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use multidimensional array approach
% Mttx = Mtt( :, 1 );
% Mtty = Mtt( :, 2 );
% Mttz = Mtt( :, 3 );
% 
% Mttp1x = Mttp1( :, 1 );
% Mttp1y = Mttp1( :, 2 );
% Mttp1z = Mttp1( :, 3 );
% 
% initSt.B( :, 1 ) = wv.b1preal * wv.breal( :, tt ) - wv.b1pimag * wv.bimag( :, tt );
% initSt.B( :, 2 ) = wv.b1pimag * wv.breal( :, tt ) + wv.b1preal * wv.bimag( :, tt );
% initSt.B( :, 3 ) = wv.bzsens * wv.Shim( :, tt ) + wv.pos * wv.Grad( :, tt ) + wv.db0;
% 
% magBtt = sqrt( initSt.B( :, 1 ).^2 + initSt.B( :, 2 ).^2 + initSt.B( :, 3 ).^2 );
%
% dt = wv.dtvec( tt );
%
% phi = - ( wv.gyro * dt ) * magBtt; 
% 
% ux = initSt.B( :, 1 ) ./ magBtt;
% uy = initSt.B( :, 2 ) ./ magBtt;
% uz = initSt.B( :, 3 ) ./ magBtt;
% 
% % Precompute values
% uxuy = ux .* uy;
% uxuz = ux .* uz;
% uyuz = uy .* uz;
% 
% uxsq = ux.^2;
% uysq = uy.^2;
% uzsq = uz.^2;
% 
% cosphi = cos( phi );
% onemcosphi = ( 1 - cosphi );
% 
% sinphi = sin( phi );
% 
% udotM = ux .* Mttx + uy .* Mtty + uz .* Mttz;
% 
% % start assembling gradients
% initSt.partialfpartialu( :, 1 ) = onemcosphi .* ( udotM + ux .* Mttx );
% initSt.partialfpartialu( :, 4 ) = onemcosphi .* ( ux .* Mtty ) + sinphi .* Mttz;
% initSt.partialfpartialu( :, 7 ) = onemcosphi .* ( ux .* Mttz ) - sinphi .* Mtty;
% initSt.partialfpartialu( :, 2 ) = onemcosphi .* ( uy .* Mttx ) - sinphi .* Mttz;
% initSt.partialfpartialu( :, 5 ) = onemcosphi .* ( udotM + uy .* Mtty );
% initSt.partialfpartialu( :, 8 ) = onemcosphi .* ( uy .* Mttz ) + sinphi .* Mttx;
% initSt.partialfpartialu( :, 3 ) = onemcosphi .* ( uz .* Mttx ) + sinphi .* Mtty;
% initSt.partialfpartialu( :, 6 ) = onemcosphi .* ( uz .* Mtty ) - sinphi .* Mttx;
% initSt.partialfpartialu( :, 9 ) = onemcosphi .* ( udotM + uz .* Mttz );
% 
% initSt.partialfpartialphi( :, 1 ) = - uz .* Mttp1y + uy .* Mttp1z;
% initSt.partialfpartialphi( :, 2 ) = + uz .* Mttp1x - ux .* Mttp1z;
% initSt.partialfpartialphi( :, 3 ) = - uy .* Mttp1x + ux .* Mttp1y;
% 
% initSt.partialupartialB( :, 1 ) = ( 1 - uxsq ) ./  magBtt;
% initSt.partialupartialB( :, 4 ) = ( - uxuy ) ./  magBtt;
% initSt.partialupartialB( :, 7 ) = ( - uxuz ) ./  magBtt;
% initSt.partialupartialB( :, 5 ) = ( 1 - uysq ) ./  magBtt;
% initSt.partialupartialB( :, 8 ) = ( - uyuz ) ./  magBtt;
% initSt.partialupartialB( :, 9 ) = ( 1 - uzsq ) ./  magBtt;
% initSt.partialupartialB( :, 2 ) = initSt.partialupartialB( :, 4 );
% initSt.partialupartialB( :, 3 ) = initSt.partialupartialB( :, 7 );
% initSt.partialupartialB( :, 6 ) = initSt.partialupartialB( :, 8 );
% 
% initSt.partialphipartialB( :, 1 ) = ( -wv.gyro * dt ) * ux;
% initSt.partialphipartialB( :, 2 ) = ( -wv.gyro * dt ) * uy;
% initSt.partialphipartialB( :, 3 ) = ( -wv.gyro * dt ) * uz;
% 
% initSt.partialfpartialB( :, 1 ) = ...
%     initSt.partialfpartialu( :, 1 ) .* initSt.partialupartialB( :, 1 ) +...
%     initSt.partialfpartialu( :, 4 ) .* initSt.partialupartialB( :, 2 ) +...
%     initSt.partialfpartialu( :, 7 ) .* initSt.partialupartialB( :, 3 ) +...
%     initSt.partialfpartialphi( :, 1 ) .* initSt.partialphipartialB( :, 1 );
% initSt.partialfpartialB( :, 4 ) = ...
%     initSt.partialfpartialu( :, 1 ) .* initSt.partialupartialB( :, 4 ) +...
%     initSt.partialfpartialu( :, 4 ) .* initSt.partialupartialB( :, 5 ) +...
%     initSt.partialfpartialu( :, 7 ) .* initSt.partialupartialB( :, 6 ) +...
%     initSt.partialfpartialphi( :, 1 ) .* initSt.partialphipartialB( :, 2 );
% initSt.partialfpartialB( :, 7 ) = ...
%     initSt.partialfpartialu( :, 1 ) .* initSt.partialupartialB( :, 7 ) +...
%     initSt.partialfpartialu( :, 4 ) .* initSt.partialupartialB( :, 8 ) +...
%     initSt.partialfpartialu( :, 7 ) .* initSt.partialupartialB( :, 9 ) +...
%     initSt.partialfpartialphi( :, 1 ) .* initSt.partialphipartialB( :, 3 );
% initSt.partialfpartialB( :, 2 ) = ...
%     initSt.partialfpartialu( :, 2 ) .* initSt.partialupartialB( :, 1 ) +...
%     initSt.partialfpartialu( :, 5 ) .* initSt.partialupartialB( :, 2 ) +...
%     initSt.partialfpartialu( :, 8 ) .* initSt.partialupartialB( :, 3 ) +...
%     initSt.partialfpartialphi( :, 2 ) .* initSt.partialphipartialB( :, 1 );
% initSt.partialfpartialB( :, 5 ) = ...
%     initSt.partialfpartialu( :, 2 ) .* initSt.partialupartialB( :, 4 ) +...
%     initSt.partialfpartialu( :, 5 ) .* initSt.partialupartialB( :, 5 ) +...
%     initSt.partialfpartialu( :, 8 ) .* initSt.partialupartialB( :, 6 ) +...
%     initSt.partialfpartialphi( :, 2 ) .* initSt.partialphipartialB( :, 2 );
% initSt.partialfpartialB( :, 8 ) = ...
%     initSt.partialfpartialu( :, 2 ) .* initSt.partialupartialB( :, 7 ) +...
%     initSt.partialfpartialu( :, 5 ) .* initSt.partialupartialB( :, 8 ) +...
%     initSt.partialfpartialu( :, 8 ) .* initSt.partialupartialB( :, 9 ) +...
%     initSt.partialfpartialphi( :, 2 ) .* initSt.partialphipartialB( :, 3 );
% initSt.partialfpartialB( :, 3 ) = ...
%     initSt.partialfpartialu( :, 3 ) .* initSt.partialupartialB( :, 1 ) +...
%     initSt.partialfpartialu( :, 6 ) .* initSt.partialupartialB( :, 2 ) +...
%     initSt.partialfpartialu( :, 9 ) .* initSt.partialupartialB( :, 3 ) +...
%     initSt.partialfpartialphi( :, 3 ) .* initSt.partialphipartialB( :, 1 );
% initSt.partialfpartialB( :, 6 ) = ...
%     initSt.partialfpartialu( :, 3 ) .* initSt.partialupartialB( :, 4 ) +...
%     initSt.partialfpartialu( :, 6 ) .* initSt.partialupartialB( :, 5 ) +...
%     initSt.partialfpartialu( :, 9 ) .* initSt.partialupartialB( :, 6 ) +...
%     initSt.partialfpartialphi( :, 3 ) .* initSt.partialphipartialB( :, 2 );
% initSt.partialfpartialB( :, 9 ) = ...
%     initSt.partialfpartialu( :, 3 ) .* initSt.partialupartialB( :, 7 ) +...
%     initSt.partialfpartialu( :, 6 ) .* initSt.partialupartialB( :, 8 ) +...
%     initSt.partialfpartialu( :, 9 ) .* initSt.partialupartialB( :, 9 ) +...
%     initSt.partialfpartialphi( :, 3 ) .* initSt.partialphipartialB( :, 3 );


