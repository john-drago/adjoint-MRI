function [ MArr, targSt ] = generateMArray( targSt, MSt )

[ targSt.X, targSt.Y, targSt.Z ] = ndgrid( targSt.x, targSt.y, targSt.z );

MStSize = size( MSt );

MArr = zeros( [ size( targSt.X ), 3 ] );

if length( MStSize ) == 2

    MArr( :, :, :, 1 ) = MSt( 1 );
    MArr( :, :, :, 2 ) = MSt( 2 );
    MArr( :, :, :, 3 ) = MSt( 3 );

elseif length( MStSize ) == 3

    MArr( :, :, :, 1 ) = MSt;

elseif length( MStSize ) == 4
    
    MArr( :, :, :, : ) = MSt;

end

end