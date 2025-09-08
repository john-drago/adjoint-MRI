function [ oc, pulse, opt ] = processMPpTxPulse( oc, pulse, opt )
% This function will process an MP-pTx pulse to determine

%% Process specific pulse morphology
switch char( pulse.type )
    case { 'or-blip-mp' }
        [ oc, pulse, opt ] = processMPpTxPulse_orblipmp( oc, pulse, opt );
    case { 'or-blip-mp-const' }
        [ oc, pulse, opt ] = processMPpTxPulse_orblipmpconst( oc, pulse, opt );
    case { 'or-blip-mp-fixed-phase' }
        [ oc, pulse, opt ] = processMPpTxPulse_orblipmp_fixedphase( oc, pulse, opt );
    case { 'or-blip-mp-const-fixed-phase' }
        [ oc, pulse, opt ] = processMPpTxPulse_orblipmpconst_fixedphase( oc, pulse, opt );
    otherwise
        error( "Unknown MP-pTx 'pulse.type'." )
end

end