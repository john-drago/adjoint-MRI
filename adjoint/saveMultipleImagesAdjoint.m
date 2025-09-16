function success = saveMultipleImagesAdjoint( savePath, strLabels, cellFigs, binVis, saveTypeFlag )
% This function will iterate across a set of images and save them to the
% savePath location
arguments
    savePath
    strLabels
    cellFigs
    binVis
    saveTypeFlag = 1 % 1 is png only, 2 is fig also and 3 is with eps
end

%% Check to make sure that there is the correct number of labels
numLabels = length( strLabels );
numImgs = length( cellFigs );

if numLabels ~= numImgs
    error( "Not the same amount of images and labels" );
end

%% Iterate over images and save
if ~isfolder( savePath )
    mkdir( savePath )
end

for nn = 1:numLabels
    savePathIter = fullfile( savePath, strLabels{ nn } );
    
    saveas( cellFigs{ nn }, strcat( savePathIter, '.png' ) );
    if saveTypeFlag > 1
        saveas( cellFigs{ nn }, strcat( savePathIter, '.fig' ) );
    end
    if saveTypeFlag > 2
        print( cellFigs{ nn }, strcat( savePathIter, '.eps' ), "-depsc" )
    end
    
    if ~binVis
        close( cellFigs{ nn } );
    end

end

success = true;

end