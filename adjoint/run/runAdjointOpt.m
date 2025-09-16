function opt = runAdjointOpt( opt, oc )
% This function will run the adjoint opt process using the information
% contain in the opt struct

%% Run Optimizations
optTic = tic;
opt.optStartTime = datetime;
if matches( oc.optType, "ipopt", 'ignorecase', true )
    
    opt = runOpt_IPOPT( opt );

elseif matches( oc.optType, "snopt", 'ignorecase', true )

    opt = runOpt_SNOPT( opt );

elseif matches( oc.optType, "nlopt", 'ignorecase', true )
    
    opt = runOpt_NLOPT( opt );

elseif matches( oc.optType, ["fmin"; "f-min" ], 'ignorecase', true )

    opt = runOpt_FMIN( opt );

elseif matches( oc.optType, ["ga"; "genetic" ], 'ignorecase', true )

    opt = runOpt_GA( opt, oc );

elseif matches( oc.optType, "ga-proxy", 'ignorecase', true )
    
    opt = runOpt_GAProxy( opt, oc );

end

%% Get opt data
optToc = toc( optTic );
opt.optTime = optToc;
opt.timeStampIden = string( datetime('now', 'format','yyMMddHHmmss'));
opt.timeStamp = string( datetime );
opt.pScOpt = opt.pScOpt( : );
opt.pOpt = opt.pScOpt .* opt.scVec;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function opt = runOpt_IPOPT( opt )

% Create output struct
output = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run Optimization
[output.pScOpt, output.output] = ipopt(...
    opt.pSc0, opt.ipopt.funcs, opt.ipopt.options );

[ ~, ~, feasible ] =...
    ipoptObjectiveBestOptWrapper( output.pScOpt, opt.ipopt.funcs.objective, opt.ipopt.funcs, opt.ipopt.options );
clear ipoptObjectiveBestOptWrapper;

if feasible.bestCost < inf
    output.output.objective = feasible.bestCost;
    output.output.bestIter = feasible.bestIter;
    output.output.bestCost = feasible.bestCost;
    output.pScOpt = feasible.pScOpt;
end

output.exitflag = output.output.status;
output.fval = output.output.objective;

opt.output = output;
opt.pScOpt = output.pScOpt( : );

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

if opt.trackConvergence
    [ ~, ~, trackConvSt ] =...
        costFnTrackConvergenceWrapper( opt.pScOpt, opt.ipopt.funcs.objective, opt, opt.constrTolSave );
    opt.fval_conv_iters = trackConvSt.fval_conv_iters;
    opt.optTime_conv_iters = trackConvSt.optTime_conv_iters;
    opt.funccount_conv_iters = trackConvSt.funccount_conv_iters;
    opt.fvalraw_conv_iters = trackConvSt.fvalraw_conv_iters;
    opt.const_conv_iters = trackConvSt.const_conv_iters;

    clear costFnTrackConvergenceWrapper;
end

if opt.trackDecisionVariables
    [ ~, ~, trackDecSt ] =...
        costFnTrackDecisionVariablesWrapper( opt.pScOpt, opt.ipopt.funcs.objective, opt, opt.constrTolSave );
    opt.fval_vars_iters = trackDecSt.fval_vars_iters;
    opt.optTime_vars_iters = trackDecSt.optTime_vars_iters;
    opt.funccount_vars_iters = trackDecSt.funccount_vars_iters;
    opt.pSc_vars_iters = trackDecSt.xtrack_vars_iters;
    opt.fvalraw_vars_iters = trackDecSt.fvalraw_vars_iters;
    opt.const_vars_iters = trackDecSt.const_vars_iters;
    
    clear costFnTrackDecisionVariablesWrapper;

    opt.p_vars_iters = opt.pSc_vars_iters .* opt.scVec;
end
if opt.trackConvergence || opt.trackDecisionVariables
    clear costFnTrackConvergenceWrapper;
    clear costFnTrackDecisionVariablesWrapper;
    clear ipoptObjectiveBestOptWrapper;
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function opt = runOpt_SNOPT( opt )

% Create output struct
output = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run Optimization
[output.pScOpt, output.fval, output.exitflag, output.output ] = snsolve(...
    opt.costFnOpt, double(opt.pSc0),...
    opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlconSNOPT,...
    opt.snopt );

if opt.trackConvergence
    [ ~, ~, trackConvSt ] =...
        costFnTrackConvergenceWrapper( output.pScOpt, opt.costFnOpt, opt, opt.constrTolSave );
    opt.fval_conv_iters = trackConvSt.fval_conv_iters;
    opt.optTime_conv_iters = trackConvSt.optTime_conv_iters;
    opt.funccount_conv_iters = trackConvSt.funccount_conv_iters;
    opt.fvalraw_conv_iters = trackConvSt.fvalraw_conv_iters;
    opt.const_conv_iters = trackConvSt.const_conv_iters;

    clear costFnTrackConvergenceWrapper;

    if trackConvSt.bestCost < inf
        output.output.objective = trackConvSt.bestCost;
        output.fval = output.output.objective;
        output.output.bestIter = trackConvSt.bestIter;
        output.output.bestCost = trackConvSt.bestCost;
        output.pScOpt = trackConvSt.pScOpt;
    end

end

opt.output = output;
opt.pScOpt = output.pScOpt(:);

% delete SPECS file if not saving
if ~opt.saveResult
    delete( opt.snopt.specsfile );
    rmdir( fileparts( opt.snopt.specsfile ) );
end

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

if opt.trackDecisionVariables
    [ ~, ~, trackDecSt ] =...
        costFnTrackDecisionVariablesWrapper( opt.pScOpt, opt.costFnOpt, opt, opt.constrTolSave );
    opt.fval_vars_iters = trackDecSt.fval_vars_iters;
    opt.optTime_vars_iters = trackDecSt.optTime_vars_iters;
    opt.funccount_vars_iters = trackDecSt.funccount_vars_iters;
    opt.pSc_vars_iters = trackDecSt.xtrack_vars_iters;
    opt.fvalraw_vars_iters = trackDecSt.fvalraw_vars_iters;
    opt.const_vars_iters = trackDecSt.const_vars_iters;
    
    clear costFnTrackDecisionVariablesWrapper;

    opt.p_vars_iters = opt.pSc_vars_iters .* opt.scVec;
end
if opt.trackConvergence || opt.trackDecisionVariables
    clear costFnTrackConvergenceWrapper;
    clear costFnTrackDecisionVariablesWrapper;
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function opt = runOpt_NLOPT( opt )

% Create output struct
output = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run Optimization
[ output.pScOpt, output.fval, output.exitflag ] = nlopt_optimize(...
    opt.nlopt, double( opt.pSc0 ) );

if opt.trackConvergence
    [ ~, ~, trackConvSt ] =...
        costFnTrackConvergenceWrapper( output.pScOpt, opt.nlopt.min_objective, opt, opt.constrTolSave );
    opt.fval_conv_iters = trackConvSt.fval_conv_iters;
    opt.optTime_conv_iters = trackConvSt.optTime_conv_iters;
    opt.funccount_conv_iters = trackConvSt.funccount_conv_iters;
    opt.fvalraw_conv_iters = trackConvSt.fvalraw_conv_iters;
    opt.const_conv_iters = trackConvSt.const_conv_iters;

    clear costFnTrackConvergenceWrapper;

    if trackConvSt.bestCost < inf
        output.output.objective = trackConvSt.bestCost;
        output.fval = output.output.objective;
        output.output.bestIter = trackConvSt.bestIter;
        output.output.bestCost = trackConvSt.bestCost;
        output.pScOpt = trackConvSt.pScOpt;
    end
end

opt.output = output;
opt.pScOpt = output.pScOpt(:);

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

if opt.trackDecisionVariables
    [ ~, ~, trackDecSt ] =...
        costFnTrackDecisionVariablesWrapper( opt.pScOpt, opt.nlopt.min_objective, opt, opt.constrTolSave );
    opt.fval_vars_iters = trackDecSt.fval_vars_iters;
    opt.optTime_vars_iters = trackDecSt.optTime_vars_iters;
    opt.funccount_vars_iters = trackDecSt.funccount_vars_iters;
    opt.pSc_vars_iters = trackDecSt.xtrack_vars_iters;
    opt.fvalraw_vars_iters = trackDecSt.fvalraw_vars_iters;
    opt.const_vars_iters = trackDecSt.const_vars_iters;
    
    clear costFnTrackDecisionVariablesWrapper;

    opt.p_vars_iters = opt.pSc_vars_iters .* opt.scVec;
end

if opt.trackConvergence || opt.trackDecisionVariables
    clear costFnTrackConvergenceWrapper;
    clear costFnTrackDecisionVariablesWrapper;
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function opt = runOpt_FMIN( opt )

% % Create cost function
% opt.costFnOpt = @( pSc ) runAdjointCostFunction( pSc( : ), opt );

% Create output struct
output = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run Optimization
[output.pScOpt, output.fval, output.exitflag, output.output ] = fmincon(...
    opt.costFnOpt, double(opt.pSc0),...
    opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon,...
    opt.fminopt );

if opt.trackConvergence
    [ ~, ~, trackConvSt ] =...
        costFnTrackConvergenceWrapper( output.pScOpt, opt.costFnOpt, opt, opt.constrTolSave );
    opt.fval_conv_iters = trackConvSt.fval_conv_iters;
    opt.optTime_conv_iters = trackConvSt.optTime_conv_iters;
    opt.funccount_conv_iters = trackConvSt.funccount_conv_iters;
    opt.fvalraw_conv_iters = trackConvSt.fvalraw_conv_iters;
    opt.const_conv_iters = trackConvSt.const_conv_iters;

    clear costFnTrackConvergenceWrapper;
    
    if isempty( output.output.bestfeasible )
        output.output.bestfeasible.fval = inf;
    end
    if trackConvSt.bestCost < output.output.bestfeasible.fval
        output.output.objective = trackConvSt.bestCost;
        output.fval = output.output.objective;
        output.output.bestIter = trackConvSt.bestIter;
        output.output.bestCost = trackConvSt.bestCost;
        output.pScOpt = trackConvSt.pScOpt;
    else
        output.fval = output.output.bestfeasible.fval;
        output.pScOpt = output.output.bestfeasible.x;
    end
end

opt.output = output;
opt.pScOpt = output.pScOpt(:);

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

if opt.trackDecisionVariables
    [ ~, ~, trackDecSt ] =...
        costFnTrackDecisionVariablesWrapper( opt.pScOpt, opt.costFnOpt, opt, opt.constrTolSave );
    opt.fval_vars_iters = trackDecSt.fval_vars_iters;
    opt.optTime_vars_iters = trackDecSt.optTime_vars_iters;
    opt.funccount_vars_iters = trackDecSt.funccount_vars_iters;
    opt.pSc_vars_iters = trackDecSt.xtrack_vars_iters;
    opt.fvalraw_vars_iters = trackDecSt.fvalraw_vars_iters;
    opt.const_vars_iters = trackDecSt.const_vars_iters;
    
    clear costFnTrackDecisionVariablesWrapper;

    opt.p_vars_iters = opt.pSc_vars_iters .* opt.scVec;
end

if opt.trackConvergence || opt.trackDecisionVariables
    clear costFnTrackConvergenceWrapper;
    clear costFnTrackDecisionVariablesWrapper;
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function opt = runOpt_GAProxy( opt, oc )

% Set seed if necessary
if isfield( oc, 'rngstate')
    rng( oc.rngstate.Seed, oc.rngstate.Type );

    stream = RandStream.getGlobalStream;
    stream.State = oc.rngstate.State;
else
    rng( "shuffle", "mt19937ar" );
end

popt = opt.popt;
popt.gaopt = opt.gaopt;

% Create output struct
proxyOutput = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run GA Optimization
[proxyOutput.pScOpt, proxyOutput.fval, proxyOutput.exitflag, proxyOutput.output, proxyOutput.population, proxyOutput.scores] = ga(...
    popt.costFnOpt, double( popt.numVars ),...
    popt.A, popt.b, popt.Aeq, popt.beq, popt.lb, popt.ub, popt.nonlcon,...
    popt.gaopt);

proxyOutput.pScOpt = proxyOutput.pScOpt(:);
opt.pOutput = proxyOutput;

% Generate pSc0 for fmincon
wvpopt = opt.generateProxyWaveform( proxyOutput.pScOpt(:) .* popt.scVec, popt );

% Get rid of popt
% opt = rmfield( opt, 'popt' );
clear popt;

opt.pSc0 = opt.generateFullFromProxy( opt, wvpopt );

% Run Fmincon Optimization
% Run Optimization
opt.costFnOpt = @( pSc ) runAdjointCostFunction( pSc( : ), opt );

[ opt, oc ] = pSc0Process( opt, oc, oc.ensureFeasibleStart );

% create cost function
opt.costFnOpt = @( pSc ) opt.runCostFunction( pSc( : ), opt );

% determine whether want to track convergence
[ oc, opt, opt.costFnOpt ] = prepareTrackConvergence( oc, opt, opt.costFnOpt );

% determine whether want to track decision variables
% [ oc, opt, opt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt, opt.costFnOpt );
[ ~, opt, opt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt, opt.costFnOpt );

[output.pScOpt, output.fval, output.exitflag, output.output, output.lambda, output.grad] = fmincon(...
    opt.costFnOpt, double(opt.pSc0),...
    opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon,...
    opt.fminopt );

opt.output = output;
opt.pScOpt = output.pScOpt(:);

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function opt = runOpt_GA( opt, oc )

% Set seed if necessary
if isfield( oc, 'rngstate')
    rng( oc.rngstate.Seed, oc.rngstate.Type );

    stream = RandStream.getGlobalStream;
    stream.State = oc.rngstate.State;
else
    rng( "shuffle", "mt19937ar" );
end

% Create cost function
opt.costFnOpt = @( pSc ) runAdjointCostFunction( pSc( : ), opt );

% Create output struct
output = struct;

% Starting opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Starting opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

% Run Optimization
[output.pScOpt, output.fval, output.exitflag, output.output, output.population, output.scores] = ga(...
    opt.costFnOpt, double( opt.numVars ),...
    opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon,...
    opt.gaopt);

opt.output = output;
opt.pScOpt = output.pScOpt(:);

% Finished opt print statement
fprintf( "\n---------------------------------------------------------\n" )
fprintf( "Finishing opt:\t%s", string(datetime) )
fprintf( "\n---------------------------------------------------------\n" )

end
% ----------------------------------------------------------------------- %