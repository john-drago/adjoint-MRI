function [ val ] = calc3DMetrics( val, Mopt, FAopt, Mtarg, FAtarg, iden, cum_numPosSubj )

magNRMSE_Mxy = zeros( val.numSubj, 1 );
NRMSE_Mxy = zeros( val.numSubj, 1 );
magNRMSE_FA = zeros( val.numSubj, 1 );
meanFA = zeros( val.numSubj, 1 );
stdFA = zeros( val.numSubj, 1 );
nMaxMin = zeros( val.numSubj, 1 );
magCoeffVar = zeros( val.numSubj, 1 );

for nn = 1:val.numSubj
    if nn == 1
        optIdxs = 1:( cum_numPosSubj(nn));
    else
        optIdxs = ( cum_numPosSubj(nn-1)+1 ):( cum_numPosSubj(nn) );
    end

    [ magNRMSE_Mxy( nn ), NRMSE_Mxy( nn ), magNRMSE_FA( nn ),...
        meanFA( nn ), stdFA( nn ), nMaxMin( nn ), magCoeffVar( nn ) ] = ...
        calculateFAMtargMetrics( Mopt( optIdxs, : ), FAopt( optIdxs ), Mtarg( optIdxs, : ), FAtarg( optIdxs ) );
end

if strcmpi( iden, "" )
    metricIden = "";
elseif strcmpi( iden, "BCHP" )
    metricIden = "BCHP_";
end

val.(strcat( metricIden, "magNRMSE_Mxy" )) = magNRMSE_Mxy;
val.(strcat( metricIden, "NRMSE_Mxy" )) = NRMSE_Mxy;
val.(strcat( metricIden, "magNRMSE_FA" )) = magNRMSE_FA;
val.(strcat( metricIden, "meanFA" )) = meanFA;
val.(strcat( metricIden, "stdFA" )) = stdFA;
val.(strcat( metricIden, "nMaxMin" )) = nMaxMin;
val.(strcat( metricIden, "magCoeffVar" )) = magCoeffVar;

end