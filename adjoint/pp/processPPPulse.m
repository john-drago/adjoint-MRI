function [ oc, pulse, opt ] = processPPPulse( oc, pulse, opt )
% This function will process an optimal control pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processPPPulse_base( oc, pulse, opt );
    otherwise
        error( "Unknown piecewise polynomial 'pulse.type'." )
end

end