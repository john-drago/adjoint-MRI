function [ meanFA, FABCHP, MBCHP ] = getMeanFlipAngle( BCmag, bSensVec, M0vec, gyro, tBCHP )
bHP = bSensVec;
bHP( :, 1 ) = bHP( :, 1 ) * BCmag;
bHP( :, 2 ) = bHP( :, 2 ) * BCmag;
bHPMag = sqrt( bHP(:,1).^2 + bHP(:,2).^2 + bHP(:,3).^2 );
bHPUnit = bHP ./ bHPMag;
bHPUnit( isnan(bHPUnit) ) = 0;
thetaHP = -gyro * tBCHP * bHPMag;
MBCHP = calcRRF_vec( M0vec, bHPUnit,  thetaHP );
FABCHP = real( acosd( MBCHP(:,3) ) );

meanFA = mean( FABCHP );
end