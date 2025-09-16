function [ A_shimAcc, b_shimAcc ] = constraintPWCPulseShimAccel( opt, shimAccelConstr )
% This function will generate the matrices used for the accel constraints.

[ A_shimAcc, b_shimAcc ] = generateAccelMatrixArray(...
    'shim', opt.numZCoils, opt, shimAccelConstr );

end
