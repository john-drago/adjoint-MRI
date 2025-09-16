function [ magNRMSEFA_vec, meanFA_vec, stdFA_vec ] = calc_magNRMSE_FA( FAarray, FAtarg, optroiArr )

numSubj = size( FAarray, 4 );
magNRMSEFA_vec = zeros( numSubj, 1 );

if nargout > 1
    meanFA_vec = zeros( numSubj, 1 );
end
if nargout > 2
    stdFA_vec = zeros( numSubj, 1 );
end


for nn = 1:numSubj
    FASubj = FAarray( :, :, :, nn );
    FAtargSubj = FAtarg( :, :, :, nn );
    optroiSubj = optroiArr( :, :, :, nn );

    FASubjVec = FASubj( optroiSubj );
    FAtargSubjVec = FAtargSubj( optroiSubj );

    magNRMSEFA_vec( nn ) = norm( FASubjVec - FAtargSubjVec ) /...
        norm( FAtargSubjVec );

    if nargout > 1
        meanFA_vec( nn ) = mean( FASubjVec );
    end
    if nargout > 2
        stdFA_vec( nn ) = std( FASubjVec );
    end
end


end