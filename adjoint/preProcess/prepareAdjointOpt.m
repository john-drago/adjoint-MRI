function [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt )
% This function will prepare the adjoint optimization process.

if ( ~isfield( oc, 'saveResult' ) )
    oc.saveResult = false;
    opt.saveResult = false;
else
    opt.saveResult = oc.saveResult;
    opt.saveDir = oc.saveDir;
end

% specify solver type and settings
if matches( oc.optType, "ipopt", 'ignorecase', true )

    if ~isfield( oc, 'ensureFeasibleStart' )
        oc.ensureFeasibleStart = false;
    end
    
    ipopt = oc.ipopt;

    % add the optimization problem information
    [ opt, ipopt ] = prepareIPOPTFunctions( opt, ipopt );

    % process the intial guess
    [ opt, oc ] = pSc0Process( opt, oc, oc.ensureFeasibleStart ); 

    % determine whether want to track convergence
    [ oc, opt, ipopt.funcs.objective ] = prepareTrackConvergence( oc, opt, ipopt.funcs.objective );

    % determine whether want to track decision variables
    [ oc, opt, ipopt.funcs.objective ] = prepareTrackDecisionVariables( oc, opt, ipopt.funcs.objective );

    opt.ipopt = ipopt;

elseif matches( oc.optType, "snopt", 'ignorecase', true )
    
    if ~isfield( oc, 'ensureFeasibleStart' )
        oc.ensureFeasibleStart = false;
    end

    % initialize the snopt struct
    snopt = oc.snopt;

    % make SPECS file
    snopt = generateSPECSfileSNOPT( oc, pulse, snopt );

    % add the optimization problem information
    [ opt, snopt ] = prepareSNOPTFunctions( opt, snopt );

    % process the intial guess
    [ opt, oc ] = pSc0Process( opt, oc, oc.ensureFeasibleStart ); 

    % determine whether want to track convergence
    [ oc, opt, opt.costFnOpt ] = prepareTrackConvergence( oc, opt, opt.costFnOpt );

    % determine whether want to track decision variables
    [ oc, opt, opt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt, opt.costFnOpt );

    opt.snopt = snopt;

elseif matches( oc.optType, "nlopt", 'ignorecase', true )

    if ~isfield( oc, 'ensureFeasibleStart' )
        oc.ensureFeasibleStart = false;
    end

    % initialize the nlopt struct
    % this should include stopping criteria 
    nlopt = oc.nlopt; 
    
    % add the optimization problem information
    [ opt, nlopt ] = prepareNLOPTFunctions( opt, nlopt );

    % reassign to the opt struct
    [ opt, oc ] = pSc0Process( opt, oc, oc.ensureFeasibleStart );
    
    % determine whether want to track convergence
    [ oc, opt, nlopt.min_objective ] = prepareTrackConvergence( oc, opt, nlopt.min_objective );

    % determine whether want to track decision variables
    [ oc, opt, nlopt.min_objective ] = prepareTrackDecisionVariables( oc, opt, nlopt.min_objective );

    opt.nlopt = nlopt;

elseif matches( oc.optType, ["fmin"; "f-min"; "fmincon" ], 'ignorecase', true )

    if ~isfield( oc, 'ensureFeasibleStart' )

        % if matches( oc.fminopt.Algorithm, "interior-point", 'ignorecase', true )
        %     oc.ensureFeasibleStart = false;
        % elseif matches( oc.fminopt.Algorithm, [ "active-set"; "sqp" ], 'ignorecase', true )
        %     oc.ensureFeasibleStart = true;
        % end

        oc.ensureFeasibleStart = true;
    end

    opt.fminopt = oc.fminopt;

    [ opt, oc ] = pSc0Process( opt, oc, oc.ensureFeasibleStart );
    
    % create cost function
    opt.costFnOpt = @( pSc ) opt.runCostFunction( pSc( : ), opt );
    
    % determine whether want to track convergence
    [ oc, opt, opt.costFnOpt ] = prepareTrackConvergence( oc, opt, opt.costFnOpt );

    % determine whether want to track decision variables
    [ oc, opt, opt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt, opt.costFnOpt );

    % check to see if using interior-point algorithm, and if so, make the
    % constraints matrices sparse
    if matches( opt.fminopt.Algorithm, "interior-point", 'ignorecase', true )
        opt.A = sparse( opt.A );
        opt.b = sparse( opt.b );
        opt.Aeq = sparse( opt.Aeq );
        opt.beq = sparse( opt.beq );
    end

elseif matches( oc.optType, ["ga"; "genetic" ], 'ignorecase', true )

    oc.gaopt.HybridFcn = {str2func(oc.fminopt.SolverName), oc.fminopt};
    
    opt.fminopt = oc.fminopt;
    opt.gaopt = oc.gaopt;

    % create cost function
    opt.costFnOpt = @( pSc ) opt.runCostFunction( pSc( : ), opt );
    
    % determine whether want to track convergence
    [ oc, opt, opt.costFnOpt ] = prepareTrackConvergence( oc, opt, opt.costFnOpt );

    % determine whether want to track decision variables
    [ oc, opt, opt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt, opt.costFnOpt );

elseif matches( oc.optType, "ga-proxy", 'ignorecase', true )

    if ~isfield( oc, 'ensureFeasibleStart' )
        oc.ensureFeasibleStart = true;
    end
    
    oc = determineParallelWorkers( oc.numWorkers, oc );

    opt.fminopt = oc.fminopt;
    opt.gaopt = oc.gaopt;

    opt.popt = opt.generateProxy( opt, oc );

    [ oc, pulse, opt.popt ] = processAdjointFunctions( oc, pulse, opt.popt );

    opt.popt.generateWaveforms = opt.generateProxyWaveform;

    % create cost function
    opt.popt.costFnOpt = @( pSc ) runAdjointCostFunction( pSc( : ), opt.popt );
    
    % determine whether want to track convergence
    [ oc, opt.popt, opt.popt.costFnOpt ] = prepareTrackConvergence( oc, opt.popt, opt.popt.costFnOpt );

    % determine whether want to track decision variables
    [ oc, opt.popt, opt.popt.costFnOpt ] = prepareTrackDecisionVariables( oc, opt.popt, opt.popt.costFnOpt );

end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %