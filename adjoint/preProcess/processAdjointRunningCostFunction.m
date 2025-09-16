function [ opt, pulse ] = processAdjointRunningCostFunction( opt, pulse )

if ~isempty( pulse.runningCostFunction )

    if matches( pulse.runningCostFunction, 'arccos-least-squares', 'IgnoreCase', true  )

        defaultRunningCostWeightingCoefficient = 1e-2;
        [ opt, pulse ] = checkRunningCostWeightingCoefficient( opt, pulse, defaultRunningCostWeightingCoefficient );
        opt.runningCostFunction = @runArccosLeastSquaresRunningCost;
    
    elseif matches( pulse.runningCostFunction, 'RF-power-pwc', 'IgnoreCase', true  )

        defaultRunningCostWeightingCoefficient = 1e-4;
        [ opt, pulse ] = checkRunningCostWeightingCoefficient( opt, pulse, defaultRunningCostWeightingCoefficient );
        opt.runningCostFunction = @runRFPowerPWCRunningCost;

    elseif matches( pulse.runningCostFunction, 'RF-power-cheb', 'IgnoreCase', true  )

        defaultRunningCostWeightingCoefficient = 1e-4;
        [ opt, pulse ] = checkRunningCostWeightingCoefficient( opt, pulse, defaultRunningCostWeightingCoefficient );
        opt.runningCostFunction = @runRFPowerChebRunningCost;

    else
        error( "Unknown 'pulse.runningCostFunction'." );
    end

else
    opt.runningCostFunction = [];
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ opt, pulse ] = checkRunningCostWeightingCoefficient( opt, pulse, runningCostWeightingCoefficient )

if nargin < 3
    runningCostWeightingCoefficient = 1e-2;
end

if isfield( pulse, 'runningCostWeightingCoefficient' )
    opt.runningCostWeightingCoefficient = pulse.runningCostWeightingCoefficient;
else
    warning("processAdjointRunningCostFunction:NotSetRunningCostWeightingCoefficient", ...
                "Didn't set 'pulse.runningCostWeightingCoefficient'. Setting to 1e-2.");
    pulse.runningCostWeightingCoefficient = runningCostWeightingCoefficient;
    opt.runningCostWeightingCoefficient = runningCostWeightingCoefficient;
end
end