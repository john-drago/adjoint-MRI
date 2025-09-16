function [ oc, pulse, opt ] = processFourierPulse( oc, pulse, opt )
% This function will process an optimal control pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processFourierPulse_base( oc, pulse, opt );
    case { "cheb"; "chebyshev"; "Cheb"; "Chebyshev" }
        pulse.type = "cheb";
        [ oc, pulse, opt ] = processFourierPulse_cheb( oc, pulse, opt );
    otherwise
        error( "Unknown optimal control 'pulse.type'." )
end

end