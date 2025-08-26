function opt = getSliceIdxs( opt )

opt.subjSlicesGloIdxs = false( size( opt.allSlicesMIdx, 1 ), opt.numSubj );
opt.subjNSlicesGloIdxs = false( size( opt.allSlicesMIdx, 1 ), opt.numSubj );
opt.subjSlicesLocIdxs = cell( opt.numSubj, 1 );
opt.subjNSlicesLocIdxs = cell( opt.numSubj, 1 );
opt.numInSliceSubjs = zeros( opt.numSubj, 1 );
opt.numOutSliceSubjs = zeros( opt.numSubj, 1 );

for nn = 1:opt.numSubj
    if nn == 1
        subjIdxs = 1 : opt.cum_numPosSubj( nn );
    else
        subjIdxs = ( opt.cum_numPosSubj( nn-1 ) + 1 ) : opt.cum_numPosSubj( nn );
    end

    opt.subjSlicesGloIdxs( subjIdxs, nn ) = opt.allSlicesMIdx( subjIdxs, 1 );
    opt.subjNSlicesGloIdxs( subjIdxs, nn ) = ~opt.subjSlicesGloIdxs( subjIdxs, nn );

    opt.subjSlicesLocIdxs{ nn, 1 } = opt.subjSlicesGloIdxs( subjIdxs, nn );
    opt.subjNSlicesLocIdxs{ nn, 1 } = opt.subjNSlicesGloIdxs( subjIdxs, nn );

    opt.numInSliceSubjs( nn, 1 ) = sum( opt.subjSlicesLocIdxs{ nn, 1 }, 1 );
    opt.numOutSliceSubjs( nn, 1 ) = size( opt.subjSlicesLocIdxs{ nn, 1 }, 1 ) - opt.numInSliceSubjs( nn, 1 );

end

end