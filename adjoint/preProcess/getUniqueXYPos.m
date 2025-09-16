function [ uniqueXYPos, XYPosToUniqueXYPos, uniqueXYPosToXYPos, numUniqueXYPos  ] = ...
    getUniqueXYPos( pos )

XYPos = pos( :, 1:2 );
[ uniqueXYPos, XYPosToUniqueXYPos, uniqueXYPosToXYPos ] =...
    unique( XYPos, 'rows' );
XYPosToUniqueXYPos = uint32( XYPosToUniqueXYPos );
uniqueXYPosToXYPos = uint32( uniqueXYPosToXYPos );
numUniqueXYPos = size( XYPosToUniqueXYPos, 1 );
% uniqueIdxs = zeros( numPos, numUniqueXYPos );
%
% for uu = uint32( 1:numUniqueXYPos )
%     uniqueIdxs( :, uu ) = ( uniqueXYPosToXYPos == uu );
% end


end