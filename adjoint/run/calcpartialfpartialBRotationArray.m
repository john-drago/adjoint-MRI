function partialfpartialB = calcpartialfpartialBRotationArray( partialfpartialB, wv, Marray )
% This function will calculate the partialfpartialB assuming the use of a
% rotation matrix for forward model integration.


% get Mnn and Mnnp1 for later gradient calculations
Marraynn = permute( Marray( :, :, 1:(end-1) ), [ 1, 3, 2 ] );

% determine where the rotating frame updates occur
if ~wv.constantRotatingFrame
    dwxy_tol = 1e-1;
    dwxydiffs = diff( wv.dwxyvec );
    dwxydiffslocm1 = find( abs( dwxydiffs ) > dwxy_tol );
    dwxydiffsloc = dwxydiffslocm1 + 1;
    
    dwxy_rot = dwxydiffs( dwxydiffslocm1 );
    t_rot = wv.tvec( dwxydiffslocm1 ) + 0.5 * wv.dtvec( dwxydiffslocm1 );
    dwxyt_rot = dwxy_rot .* t_rot;

    cosdwxyt_rot = repmat( reshape( cos( dwxyt_rot ), [ 1, 1, length(t_rot) ] ), [ wv.numPos, 1 ] );
    sindwxyt_rot = repmat( reshape( sin( dwxyt_rot ), [ 1, 1, length(t_rot) ] ), [ wv.numPos, 1 ] );

    Marray( :, :, dwxydiffsloc ) = cat( 2, ...
        Marray( :, 1, dwxydiffsloc ) .* cosdwxyt_rot - Marray( :, 2, dwxydiffsloc ) .* sindwxyt_rot,...
        Marray( :, 2, dwxydiffsloc ) .* cosdwxyt_rot + Marray( :, 1, dwxydiffsloc ) .* sindwxyt_rot,...
        Marray( :, 3, dwxydiffsloc ) );

end

Marraynnp1 = permute( Marray( :, :, 2:(end) ), [ 1, 3, 2 ] );

clear Marray;

Bx = wv.b1preal * wv.breal - wv.b1pimag * wv.bimag;
By = wv.b1pimag * wv.breal + wv.b1preal * wv.bimag;

if ( wv.numZCoils > 0 )
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

magBnn = sqrt( Bx.^2 + By.^2 + Bz.^2 );

phi = - magBnn .* ( wv.gyro * wv.dtvec );

ux = Bx ./ magBnn;
uy = By ./ magBnn;
uz = Bz ./ magBnn;

% Precompute values
% uxuy = ux .* uy;
% uxuz = ux .* uz;
% uyuz = uy .* uz;

% uxsq = ux.^2;
% uysq = uy.^2;
% uzsq = uz.^2;
% cosphi = cos( phi );
% onemcosphi = ( 1 - cosphi );

onemcosphi = ( 1 - cos( phi ) );
sinphi = sin( phi );

clear Bx By Bz phi;

Marraynnx = Marraynn( :, :, 1 );
Marraynny = Marraynn( :, :, 2 );
Marraynnz = Marraynn( :, :, 3 );

clear Marraynn;

udotM = ux .* Marraynnx + uy .* Marraynny + uz .* Marraynnz;

% start assembling gradients
partialfpartialu_1 = onemcosphi .* ( udotM + ux .* Marraynnx );
partialfpartialu_4 = onemcosphi .* ( ux .* Marraynny ) + sinphi .* Marraynnz;
partialfpartialu_7 = onemcosphi .* ( ux .* Marraynnz ) - sinphi .* Marraynny;
partialfpartialu_2 = onemcosphi .* ( uy .* Marraynnx ) - sinphi .* Marraynnz;
partialfpartialu_5 = onemcosphi .* ( udotM + uy .* Marraynny );
partialfpartialu_8 = onemcosphi .* ( uy .* Marraynnz ) + sinphi .* Marraynnx;
partialfpartialu_3 = onemcosphi .* ( uz .* Marraynnx ) + sinphi .* Marraynny;
partialfpartialu_6 = onemcosphi .* ( uz .* Marraynny ) - sinphi .* Marraynnx;
partialfpartialu_9 = onemcosphi .* ( udotM + uz .* Marraynnz );

clear onemcosphi sinphi udotM Marraynnx Marraynny Marraynnz;

partialfpartialphi_1 = - uz .* Marraynnp1( :, :, 2 ) + uy .* Marraynnp1( :, :, 3 );
partialfpartialphi_2 = + uz .* Marraynnp1( :, :, 1 ) - ux .* Marraynnp1( :, :, 3 );
partialfpartialphi_3 = - uy .* Marraynnp1( :, :, 1 ) + ux .* Marraynnp1( :, :, 2 );

clear Marraynnp1;

partialupartialB_1 = ( 1 - ( ux.^2 ) ) ./  magBnn;
partialupartialB_4 = ( - ( ux .* uy ) ) ./  magBnn;
partialupartialB_7 = ( - ( ux .* uz ) ) ./  magBnn;
partialupartialB_5 = ( 1 - ( uy.^2 ) ) ./  magBnn;
partialupartialB_8 = ( - ( uy .* uz ) ) ./  magBnn;
partialupartialB_9 = ( 1 - ( uz.^2 ) ) ./  magBnn;
% partialupartialB_2 = partialupartialB_4;
% partialupartialB_3 = partialupartialB_7;
% partialupartialB_6 = partialupartialB_8;

clear magBnn;

partialphipartialB_1 = ux .* ( -wv.gyro * wv.dtvec );
partialphipartialB_2 = uy .* ( -wv.gyro * wv.dtvec );
partialphipartialB_3 = uz .* ( -wv.gyro * wv.dtvec );

partialfpartialB( :, :, 1 ) = ...
    partialfpartialu_1 .* partialupartialB_1 +...
    partialfpartialu_4 .* partialupartialB_4 +... % should be partialfpartialu_4 .* partialupartialB_2 +...
    partialfpartialu_7 .* partialupartialB_7 +... % should be partialfpartialu_7 .* partialupartialB_3 +...
    partialfpartialphi_1 .* partialphipartialB_1;


partialfpartialB( :, :, 4 ) = ...
    partialfpartialu_1 .* partialupartialB_4 +... 
    partialfpartialu_4 .* partialupartialB_5 +...
    partialfpartialu_7 .* partialupartialB_8 +... % should be partialfpartialu_7 .* partialupartialB_6 +...
    partialfpartialphi_1 .* partialphipartialB_2;


partialfpartialB( :, :, 7 ) = ...
    partialfpartialu_1 .* partialupartialB_7 +...
    partialfpartialu_4 .* partialupartialB_8 +...
    partialfpartialu_7 .* partialupartialB_9 +...
    partialfpartialphi_1 .* partialphipartialB_3;


partialfpartialB( :, :, 2 ) = ...
    partialfpartialu_2 .* partialupartialB_1 +...
    partialfpartialu_5 .* partialupartialB_4 +... % should be partialfpartialu_5 .* partialupartialB_2 +...
    partialfpartialu_8 .* partialupartialB_7 +... % should be partialfpartialu_8 .* partialupartialB_3 +...
    partialfpartialphi_2 .* partialphipartialB_1;


partialfpartialB( :, :, 5 ) = ...
    partialfpartialu_2 .* partialupartialB_4 +...
    partialfpartialu_5 .* partialupartialB_5 +...
    partialfpartialu_8 .* partialupartialB_8 +... % should be partialfpartialu_8 .* partialupartialB_6 +...
    partialfpartialphi_2 .* partialphipartialB_2;


partialfpartialB( :, :, 8 ) = ...
    partialfpartialu_2 .* partialupartialB_7 +...
    partialfpartialu_5 .* partialupartialB_8 +...
    partialfpartialu_8 .* partialupartialB_9 +...
    partialfpartialphi_2 .* partialphipartialB_3;


partialfpartialB( :, :, 3 ) = ...
    partialfpartialu_3 .* partialupartialB_1 +...
    partialfpartialu_6 .* partialupartialB_4 +... % should be partialfpartialu_6 .* partialupartialB_2 +...
    partialfpartialu_9 .* partialupartialB_7 +... % should be partialfpartialu_9 .* partialupartialB_3 +...
    partialfpartialphi_3 .* partialphipartialB_1;


partialfpartialB( :, :, 6 ) = ...
    partialfpartialu_3 .* partialupartialB_4 +...
    partialfpartialu_6 .* partialupartialB_5 +...
    partialfpartialu_9 .* partialupartialB_8 +... % should be partialfpartialu_9 .* partialupartialB_6 +...
    partialfpartialphi_3 .* partialphipartialB_2;


partialfpartialB( :, :, 9 ) = ...
    partialfpartialu_3 .* partialupartialB_7 +...
    partialfpartialu_6 .* partialupartialB_8 +...
    partialfpartialu_9 .* partialupartialB_9 +...
    partialfpartialphi_3 .* partialphipartialB_3;

end



