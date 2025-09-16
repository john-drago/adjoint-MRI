function sliceIdx = sliceIdxNumImgs( numImgsDes, roiFOV )
numSlices = size( roiFOV,  3 );
if numSlices > numImgsDes
    sliceSpacing = numSlices / numImgsDes;
    sliceIdx = round( sliceSpacing/2 : sliceSpacing : numSlices );
else
    sliceIdx = 1:numSlices;
end
end