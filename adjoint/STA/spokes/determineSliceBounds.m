function sliceBounds = determineSliceBounds( sliceLocation, sliceThickness )

numSlices = size( sliceLocation, 1 );
sliceBounds = zeros( numSlices, 2 );
for ss = 1:numSlices
    sliceBounds( ss, : ) = sliceLocation( ss, 1 ) + [ -sliceThickness/2, sliceThickness/2 ];
end

end