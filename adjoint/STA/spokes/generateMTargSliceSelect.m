function targSt = generateMTargSliceSelect( targSt, pulse )

numDigAfterDecimal = 8;
dttol = min( [...
    min( abs( diff( targSt.x ) ) ),...
    min( abs( diff( targSt.y ) ) ),...
    min( abs( diff( targSt.z ) ) ),...
    ] ) / 10;

x = round( targSt.x, numDigAfterDecimal );
y = round( targSt.y, numDigAfterDecimal );
z = round( targSt.z, numDigAfterDecimal );

[ X_Mtarg, Y_Mtarg, Z_Mtarg ] = ndgrid( x, y, z );

MtargValsX = zeros( size( X_Mtarg ) );
MtargValsY = zeros( size( Y_Mtarg ) );
MtargValsZ = ones( size( Z_Mtarg ) );

numSlices = size( pulse.sliceBounds, 1 );

for nn = 1:numSlices

    MtargValsX( ( Z_Mtarg >= ( pulse.sliceBounds( nn, 1 ) - dttol ) ) & ...
        ( Z_Mtarg <= ( pulse.sliceBounds( nn, 2 ) + dttol ) ) ) = sind( pulse.targFAVal );
    MtargValsZ( ( Z_Mtarg >= ( pulse.sliceBounds( nn, 1 ) - dttol ) ) & ...
        ( Z_Mtarg <= ( pulse.sliceBounds( nn, 2 ) + dttol ) ) ) = cosd( pulse.targFAVal );

end

targSt.MtargVals = cat( 4, MtargValsX, MtargValsY, MtargValsZ );

end