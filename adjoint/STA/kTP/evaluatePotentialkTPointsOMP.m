function [ K, AE, bE, RF_best, resid_best, idx_best ] =...
    evaluatePotentialkTPointsOMP( AE, bE, RF0p, K, kTTime, kT_trial, opt, pulse, fminopt )

if nargin < 9
    % initialize optimization parameters for solve
    fminopt = optimoptions( "fmincon" );
    fminopt.ConstraintTolerance = 1e-9;
    fminopt.Display = "off";
    fminopt.SpecifyConstraintGradient = true;
    fminopt.SpecifyObjectiveGradient = true;
    fminopt.StepTolerance = 1e-14;
    fminopt.FunctionTolerance = 1e-8;
    fminopt.MaxFunctionEvaluations = inf;

    fminopt.Algorithm = "active-set";
    fminopt.RelLineSrchBnd = 1e-1;
    fminopt.RelLineSrchBndDuration = inf;
    fminopt.TolConSQP = 1e-10;
    fminopt.MaxIterations = 500;
end

% get file parameters
curr_num_kTP = size( K, 2 ) + 1;

numTrials = size( kT_trial, 2 );
residNorm_trials = zeros( numTrials, 1 );

residNorm_best = inf;
idx_best = 0;

for nn = 1:numTrials
    Ann = makekTPSTAColumns(...
        complex( opt.b1preal, opt.b1pimag ), kT_trial( :, nn ), kTTime, opt, pulse );  
    AEp = [ Ann, AE ];

    [ RF, resid ] = solveOMPMatriceskTP( AEp, bE, RF0p, opt, curr_num_kTP, fminopt );

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
K = [ kT_trial( :, idx_best ), K ];

end