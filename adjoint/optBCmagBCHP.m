function [BCmagOpt, fval] = optBCmagBCHP( initBCmag, targFA_value, bSensVec, M0Vec, gyro, tBCHP )
    costFn = @(BCmag) getFlipAngleError( BCmag, targFA_value, bSensVec, M0Vec, gyro, tBCHP );
    [BCmagOpt, fval] = fminsearch( costFn, initBCmag );
end