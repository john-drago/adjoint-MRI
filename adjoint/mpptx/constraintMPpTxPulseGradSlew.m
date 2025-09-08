function [ A_gradSlew, b_gradSlew ] = constraintMPpTxPulseGradSlew( opt, gradSlewConstr )
% This function will generate the matrices used for the slew constraints
% during the RF pulses for the MPSP and ORSP.

[ A_gradSlew_x_ORSP, b_gradSlew_x_ORSP ] = generateSlewRateMatrix(...
    'Gx-ORSP', opt, gradSlewConstr );
[ A_gradSlew_y_ORSP, b_gradSlew_y_ORSP ] = generateSlewRateMatrix(...
    'Gy-ORSP', opt, gradSlewConstr );
[ A_gradSlew_z_ORSP, b_gradSlew_z_ORSP ] = generateSlewRateMatrix(...
    'Gz-ORSP', opt, gradSlewConstr );

A_gradSlew = [ A_gradSlew_x_ORSP; A_gradSlew_y_ORSP; A_gradSlew_z_ORSP ];

b_gradSlew = [ b_gradSlew_x_ORSP; b_gradSlew_y_ORSP; b_gradSlew_z_ORSP ];

end
