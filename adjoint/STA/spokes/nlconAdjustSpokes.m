function [ nonlconIneqFuncsSpokes, nonlconEqFuncsSpokes ] = nlconAdjustSpokes( opt )
% This file will add in nonlinear constraints for the spokes pulses to
% ensure compliance

%% Initialize Constraint List
nonlconIneqFuncsSpokes = {};
nonlconIneqFuncsSpokesIdx = 0;

nonlconEqFuncsSpokes = {};
nonlconEqFuncsSpokesIdx = 0;

%% Iterate over inequality constraints first
if ~isempty( opt.nlconIneqFuncs )
    for nn = 1:length( opt.nlconIneqFuncs )
        spokeFunc = opt.nlconIneqFuncs{ nn, 1 };
        spokeConstraintFunc = ...
            determineSpokeConstraintFunction( spokeFunc ); 
        if ~isempty( spokeConstraintFunc )
            nonlconIneqFuncsSpokesIdx = nonlconIneqFuncsSpokesIdx + 1;
            nonlconIneqFuncsSpokes{ nonlconIneqFuncsSpokesIdx, 1 } =...
                spokeConstraintFunc;%#ok
        end
    end
end

%% Iterate over equality constraints
if ~isempty( opt.nlconEqFuncs )
    for nn = 1:length( opt.nlconEqFuncs )
        spokeFunc = opt.nlconEqFuncs{ nn, 1 };
        spokeConstraintFunc = ...
            determineSpokeConstraintFunction( spokeFunc ); 
        if ~isempty( spokeConstraintFunc )
            nonlconEqFuncsSpokesIdx = nonlconEqFuncsSpokesIdx + 1;
            nonlconEqFuncsSpokes{ nonlconEqFuncsSpokesIdx, 1 } =...
                spokeConstraintFunc;%#ok
        end
    end
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function spokeConstraintFunc = determineSpokeConstraintFunction( spokeFunc )

spokeFuncStr = func2str( spokeFunc );

if contains( spokeFuncStr, "TotalRFPower", 'ignorecase', true )
    
    spokeConstraintFunc = @constraintSpokesTotalRFPower;

elseif contains( spokeFuncStr, "MaxRFPower", 'ignorecase', true )

    spokeConstraintFunc = @constraintSpokesMaxRFPower;

elseif contains( spokeFuncStr, "PeakLocalSAR", 'ignorecase', true )

    spokeConstraintFunc = @constraintSpokesPeakLocalSAR;

elseif contains( spokeFuncStr, "PeakGlobalSAR", 'ignorecase', true )

    spokeConstraintFunc = @constraintSpokesPeakGlobalSAR;

elseif contains( spokeFuncStr, "AvgLocalSAR", 'ignorecase', true )

    spokeConstraintFunc = @constraintSpokesAvgLocalSAR;

elseif contains( spokeFuncStr, "AvgGlobalSAR", 'ignorecase', true )

    spokeConstraintFunc = @constraintSpokesAvgGlobalSAR;

else

    spokeConstraintFunc = [];

end

end
% ----------------------------------------------------------------------- %