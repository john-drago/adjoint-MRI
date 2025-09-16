function [ A_RFslew, b_RFslew ] = constraintPWCPulseRFSlew( opt, RFSlewConstr )
% This function will generate the matrices used for the slew constraints
% during the RF pulses.

[ A_RFslew_breal, b_RFslew_breal ] = generateSlewRateMatrixArray(...
    'breal', opt.numXYCoils, opt, RFSlewConstr );
[ A_RFslew_bimag, b_RFslew_bimag ] = generateSlewRateMatrixArray(...
    'bimag', opt.numXYCoils, opt, RFSlewConstr );

A_RFslew = [ A_RFslew_breal; A_RFslew_bimag ];

b_RFslew = [ b_RFslew_breal; b_RFslew_bimag ];

end
