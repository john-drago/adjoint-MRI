function [ A_gradAcc, b_gradAcc ] = constraintMPpTxPulseGradAccel( opt, gradAccelConstr )
% This function will generate the matrices used for the slew constraints
% during the RF pulses for the MPSP and ORSP.

[ A_gradAcc_x_ORSP, b_gradAcc_x_ORSP ] = generateAccelMatrix(...
    'Gx-ORSP', opt, gradAccelConstr );
[ A_gradAcc_y_ORSP, b_gradAcc_y_ORSP ] = generateAccelMatrix(...
    'Gy-ORSP', opt, gradAccelConstr );
[ A_gradAcc_z_ORSP, b_gradAcc_z_ORSP ] = generateAccelMatrix(...
    'Gz-ORSP', opt, gradAccelConstr );

A_gradAcc = [ A_gradAcc_x_ORSP; A_gradAcc_y_ORSP; A_gradAcc_z_ORSP ];

b_gradAcc = [ b_gradAcc_x_ORSP; b_gradAcc_y_ORSP; b_gradAcc_z_ORSP ];

end
