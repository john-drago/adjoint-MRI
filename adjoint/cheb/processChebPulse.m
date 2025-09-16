function [ oc, pulse, opt ] = processChebPulse( oc, pulse, opt )
% This function will process an optimal control pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processChebPulse_base( oc, pulse, opt );
    otherwise
        error( "Unknown optimal control 'pulse.type'." )
end

end