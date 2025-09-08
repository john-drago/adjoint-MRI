function [ oc, pulse, opt ] = processSPINSPulse( oc, pulse, opt )
% This function will process a SPINS to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processSPINSPulse_base( oc, pulse, opt );
    otherwise
        error( "Unknown SPINS 'pulse.type'." )
end

end