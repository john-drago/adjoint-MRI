function [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt )
% This function will assign the correct cost functions to the opt struct
% based on definitions.

%% Determine whether to update magnetization target
if opt.convertMtargAwayLarmor && ( opt.dwxy ~= 0 )
    opt.Mtarg = convertMAwayLarmor( opt.Mtarg, opt, opt.pulseLength );
end

%% Determine type of pulse
if isfield( pulse, 'excitationType' )
    if matches( pulse.excitationType, [ "slice-selective"; "sliceselective" ],...
            'IgnoreCase', true )

        pulse.sliceSelective = true;
        pulse.nonSelective = false;
        [ opt, pulse ] = processSliceSelectivePulse( opt, pulse );
        pulse.numSlices = length( pulse.sliceLocation );

    elseif matches( pulse.excitationType, [ "non-selective"; "nonselective" ],...
            'IgnoreCase', true  )
        pulse.sliceSelective = false;
        pulse.nonSelective = true;
    else
        error( "Unknown 'pulse.excitationType'." );
    end
else
    pulse.sliceSelective = false;
    pulse.nonSelective = true;
end

opt.sliceSelective = pulse.sliceSelective;
opt.nonSelective = pulse.nonSelective;

%% Determine whether to perform BCHP comparison
if opt.nonSelective
    if isfield( oc, 'performBCHPcomp' )
        opt.performBCHPcomp = oc.performBCHPcomp;
    else
        opt.performBCHPcomp = false;
    end
else
    opt.performBCHPcomp = false;
end

%% Organize nonlcon constraint functions
if ( ~isempty( opt.nlconIneqFuncs ) ) || ( ~isempty( opt.nlconEqFuncs ) )
    opt.nonlcon = @( pSc ) nonlconCombine( pSc( : ), opt );
else
    opt.nonlcon = [];
end

%% Assign terminal cost functions
[ opt, pulse ] = processAdjointTerminalCostFunction( opt, pulse );

%% Assign running cost functions
[ opt, pulse ] = processAdjointRunningCostFunction( opt, pulse );

%% Assign cost function
opt.runCostFunction = @runAdjointCostFunction;

end