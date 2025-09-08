function [ oc, pulse, opt ] = processPulseName( oc, pulse, opt )
% This function will process the pulse struct to determine different
% aspects about the pulse that will be important for the optimization
% process.

%% Determine pulse type
switch char( pulse.name )
    case {"mpptx"; "mp-ptx"}
        pulse.name = "mpptx";
        [ oc, pulse, opt ] = processMPpTxPulse( oc, pulse, opt );
    case { "ktp"; "kTp"; "kt"; "kT"; "kTP" }
        pulse.name = "kTP";
        [ oc, pulse, opt ] = processkTPPulse( oc, pulse, opt );
    case { "variable-ktp"; "var-ktp"; "varktp"; "variablektp"; "var-kTP"; "variable-kTP"; "varkTP"; "variablekTP" }
        pulse.name = "varkTP";
        [ oc, pulse, opt ] = processvarkTPPulse( oc, pulse, opt );
    case { "PWC"; "pwc" }
        pulse.name = "pwc";
        [ oc, pulse, opt ] = processPWCPulse( oc, pulse, opt );
    case { "SPINS"; "spins" }
        pulse.name = "SPINS";
        [ oc, pulse, opt ] = processSPINSPulse( oc, pulse, opt );
    case { "cheb"; "chebyshev"; "Cheb"; "Chebyshev" }
        pulse.name = "cheb";
        [ oc, pulse, opt ] = processChebPulse( oc, pulse, opt );
    case { "four"; "fourier"; "Four"; "Fourier" }
        pulse.name = "fourier";
        [ oc, pulse, opt ] = processFourierPulse( oc, pulse, opt );
    case { "pp"; "PP" }
        pulse.name = "pp";
        [ oc, pulse, opt ] = processPPPulse( oc, pulse, opt );
    otherwise
        error("Unknown 'pulse.name'.")
end

end