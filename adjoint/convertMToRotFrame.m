function Mnn = convertMToRotFrame( MnnLarmor, dwxy, tnn )

if isgpuarray( MnnLarmor )
    dwxy = gpuArray( dwxy );
    tnn = gpuArray( tnn );
end

dwxytnn = dwxy * tnn;
cosdwxytnn = cos( dwxytnn );
sindwxytnn = sin( dwxytnn );

Mnn = cat( 2, ...
    MnnLarmor( :, 1, : ) * cosdwxytnn + MnnLarmor( :, 2, : ) * sindwxytnn,...
    MnnLarmor( :, 2, : ) * cosdwxytnn - MnnLarmor( :, 1, : ) * sindwxytnn,...
    MnnLarmor( :, 3, : )...
    );

end