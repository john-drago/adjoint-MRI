function [ val, opt ] = postProcessAdjointShim( val, opt, pulse, fields )
% This file with process the optimization results from the adjoint method.

[ val, opt ] = postProcessAdjointBase( val, opt, pulse, fields );

end


