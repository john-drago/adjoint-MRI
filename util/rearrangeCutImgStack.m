function img = rearrangeCutImgStack(imgStack, percHorzReduce, percVertReduce, numRows, numCols )
arguments
    imgStack
    percHorzReduce = 0.25;
    percVertReduce = 0.05;
    numRows = []
    numCols = []
end
% Assumes an image stack with three dimensions. The first two dimensions
% are the x and y dimensions and the third dimension is the z dimension.
%
% Will plot z slices in a mosaic

sliceSize = size(imgStack, [1 2]);
horzIdxs = [1+floor(sliceSize(1)*percHorzReduce),...
    sliceSize(1)-floor(sliceSize(1)*percHorzReduce)];
horzIdxLen = horzIdxs(2)-horzIdxs(1) + 1;
vertIdxs = [1+floor(sliceSize(2)*percVertReduce),...
    sliceSize(2)-floor(sliceSize(2)*percVertReduce)];
vertIdxLen = vertIdxs(2)-vertIdxs(1) + 1;

numSlices = size(imgStack, 3);

if isempty( numRows )
    numRows = floor(sqrt(numSlices));
end
if isempty( numCols )
    numCols = ceil(numSlices/numRows);
end

numLastRow = numSlices - (numRows-1)*numCols;
totalHorzCells = numCols*horzIdxLen;
totalCellsLastRow = numLastRow*horzIdxLen;
leftOverLastRow = totalHorzCells - totalCellsLastRow;
lastRowStartingPoint = floor(leftOverLastRow/2);

img = zeros(numCols*horzIdxLen, numRows*vertIdxLen);

imgCtr = 0;
for ss = numSlices:-1:1
    imgCtr = imgCtr + 1;

    rowIdx = ceil( imgCtr/numCols );
    colIdx = imgCtr - numCols*(rowIdx-1);

%     [imgCtr, rowIdx, colIdx]

    if rowIdx == numRows
        imgRowsIdx = lastRowStartingPoint + (colIdx-1)*horzIdxLen + (1:horzIdxLen);
        imgColsIdx = numRows*vertIdxLen - rowIdx*vertIdxLen + (1:vertIdxLen);
    else
        imgColsIdx = numRows*vertIdxLen - rowIdx*vertIdxLen + (1:vertIdxLen);
        imgRowsIdx = (colIdx-1)*horzIdxLen + (1:horzIdxLen);
    end
    
    img( imgRowsIdx, imgColsIdx ) =...
         imgStack(horzIdxs(1):horzIdxs(2), vertIdxs(1):vertIdxs(2),ss);
end

end