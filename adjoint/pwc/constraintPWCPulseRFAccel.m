function [ A_RFaccel, b_RFaccel ] = constraintPWCPulseRFAccel( opt, RFAccelConstr )
% This function will generate the matrices used for the slew constraints
% during the RF pulses.

[ A_RFaccel_breal, b_RFaccel_breal ] = generateAccelMatrixArray(...
    'breal', opt.numXYCoils, opt, RFAccelConstr );
[ A_RFaccel_bimag, b_RFaccel_bimag ] = generateAccelMatrixArray(...
    'bimag', opt.numXYCoils, opt, RFAccelConstr );

A_RFaccel = [ A_RFaccel_breal; A_RFaccel_bimag ];

b_RFaccel = [ b_RFaccel_breal; b_RFaccel_bimag ];

end
