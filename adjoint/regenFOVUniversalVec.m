function [ scArr, optroiArr, valroiArr ] = regenFOVUniversalVec(...
    FOV, scalarvalvec, idxvec, optroivec, cum_numPosSubj )

numSubj = length( cum_numPosSubj );
posStart = [ 1; (cum_numPosSubj(1:(end-1))-1) ];
posEnd = cum_numPosSubj;

xlength = length(FOV.x);
ylength = length(FOV.y);
zlength = length(FOV.z);
subjArrSize = [ xlength, ylength, zlength ];

% [ I, J, K ] = ndgrid( 1:xlength, 1:ylength, 1:zlength );

scArr = zeros( xlength, ylength, zlength, numSubj );
optroiArr = false( xlength, ylength, zlength, numSubj );
valroiArr = false( xlength, ylength, zlength, numSubj );

for nn = 1:numSubj
    subjIdx = posStart(nn):posEnd(nn);

    subjValIdx = idxvec( subjIdx );
    subjOptIdx = subjValIdx( logical( optroivec( subjIdx ) ) );
    scValIdx = scalarvalvec( subjIdx );
    % [ isubjval, jsubjval, ksubjval ] = ind2sub( subjArrSize, subjValIdx );
    % [ isubjopt, jsubjopt, ksubjopt ] = ind2sub( subjArrSize, subjOptIdx );
    
    optroiSubj = false( subjArrSize );
    valroiSubj = false( subjArrSize );
    scArrSubj = zeros( subjArrSize );

    optroiSubj( subjOptIdx ) = true;
    valroiSubj( subjValIdx ) = true;

    scArrSubj( subjValIdx ) = scValIdx;
    
    optroiArr( :, :, :, nn ) = optroiSubj;
    valroiArr( :, :, :, nn ) = valroiSubj;
    scArr( :, :, :, nn ) = scArrSubj;
end

end