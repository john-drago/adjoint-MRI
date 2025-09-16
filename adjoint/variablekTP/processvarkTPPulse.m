function [ oc, pulse, opt ] = processvarkTPPulse( oc, pulse, opt )
% This function will process an MP-pTx pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processvarkTPPulse_base( oc, pulse, opt );
    otherwise
        error( "Unknown variable kT-point 'pulse.type'." )
end

end