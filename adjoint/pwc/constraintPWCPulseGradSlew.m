function [ A_gradSlew, b_gradSlew ] = constraintPWCPulseGradSlew( opt, gradSlewConstr )
% This function will generate the matrices used for the slew constraints for the gradient.

[ A_gradSlew, b_gradSlew ] = generateSlewRateMatrixArray(...
    'grad', 3, opt, gradSlewConstr );

end
