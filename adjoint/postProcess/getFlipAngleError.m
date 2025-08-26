function [ cost ] = getFlipAngleError( BCmag, targFA, bSensVec, M0Vec, gyro, tBCHP )

[ ~, FABCHP, ~ ] = getMeanFlipAngle( BCmag, bSensVec, M0Vec, gyro, tBCHP );

cost = norm( FABCHP -  targFA ) / norm(targFA);

end