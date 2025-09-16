function MnnLarmor = convertMFromRotFrame( Mnn, dwxy, tnn )

if isgpuarray( Mnn )
    dwxy = gpuArray( dwxy );
    tnn = gpuArray( tnn );
end

dwxytnn = dwxy * tnn;
cosdwxytnn = cos( dwxytnn );
sindwxytnn = sin( dwxytnn );

MnnLarmor = cat( 2, ...
    Mnn( :, 1, : ) * cosdwxytnn - Mnn( :, 2, : ) * sindwxytnn,...
    Mnn( :, 2, : ) * cosdwxytnn + Mnn( :, 1, : ) * sindwxytnn,...
    Mnn( :, 3, : )...
    );

end