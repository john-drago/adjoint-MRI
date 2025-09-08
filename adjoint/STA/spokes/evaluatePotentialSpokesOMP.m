function [ K, AE, bE, RF_best, resid_best, idx_best ] =...
    evaluatePotentialSpokesOMP( AE, bE, RF0p, K, spokeIdx, spoke_trial, opt, spokes, fminopt )

if nargin < 9
    % initialize optimization parameters for solve
    fminopt = optimoptions( "fmincon" );
    fminopt.Display = "off";
    fminopt.SpecifyConstraintGradient = true;
    fminopt.SpecifyObjectiveGradient = true;
    fminopt.StepTolerance = 1e-5;
    fminopt.FunctionTolerance = 1e-6;
    fminopt.ConstraintTolerance = 1e-14;
    fminopt.OptimalityTolerance = 1e-7;
    fminopt.MaxFunctionEvaluations = inf;

    fminopt.MaxIterations = 500;

    fminopt.Algorithm = "sqp";

    % fminopt.Algorithm = "active-set";
    % fminopt.RelLineSrchBnd = 1e-1;
    % fminopt.RelLineSrchBndDuration = inf;
    % fminopt.TolConSQP = 1e-10;
    
end

% get file parameters
curr_numSpokes = size( K, 2 ) + 1;

numTrials = size( spoke_trial, 2 );
residNorm_trials = zeros( numTrials, 1 );

residNorm_best = inf;
idx_best = 0;

for nn = 1:numTrials
    Ann = makeSpokesSTAColumns(...
        complex( opt.b1preal, opt.b1pimag ), spoke_trial( :, nn ), spokeIdx, opt, spokes );
    AEp = [ Ann, AE ];

    [ RF, resid ] = solveOMPMatricesSpokes( AEp, bE, RF0p, spokes, opt, curr_numSpokes, fminopt );

    residNorm_trials( nn ) = norm( resid, 2 );

    if residNorm_trials( nn ) < residNorm_best
        residNorm_best = residNorm_trials( nn );
        resid_best = resid;
        Ann_best = Ann;
        RF_best = RF;
        idx_best = nn;
    end

end

% Append the new columns
AE = [ Ann_best, AE ];

% Update phase
bE = abs( bE ) .* exp( 1j * angle( AE*RF_best ) );

% Update K
K = [ spoke_trial( :, idx_best ), K ];

end