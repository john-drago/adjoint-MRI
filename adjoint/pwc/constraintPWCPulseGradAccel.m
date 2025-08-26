function [ A_gradAcc, b_gradAcc ] = constraintPWCPulseGradAccel( opt, gradAccelConstr )
% This function will generate the matrices used for the slew constraints
% during the RF pulses for the MPSP and ORSP.

[ A_gradAcc, b_gradAcc ] = generateAccelMatrixArray(...
    'grad', 3, opt, gradAccelConstr );

end
