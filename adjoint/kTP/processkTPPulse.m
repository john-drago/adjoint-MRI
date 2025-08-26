function [ oc, pulse, opt ] = processkTPPulse( oc, pulse, opt )
% This function will process an MP-pTx pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'base' }
        [ oc, pulse, opt ] = processkTPPulse_base( oc, pulse, opt );
    otherwise
        error( "Unknown kT-point 'pulse.type'." )
end

end