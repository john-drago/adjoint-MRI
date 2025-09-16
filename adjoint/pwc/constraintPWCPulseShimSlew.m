function [ A_shimSlew, b_shimSlew ] = constraintPWCPulseShimSlew( opt, shimSlewConstr )
% This function will generate the matrices used for the slew constraints for the shim array.

[ A_shimSlew, b_shimSlew ] = generateSlewRateMatrixArray(...
    'shim', opt.numZCoils, opt, shimSlewConstr );

end
