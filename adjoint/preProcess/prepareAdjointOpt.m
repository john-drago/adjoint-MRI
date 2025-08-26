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

elseif matches( oc.optType, "ga-proxy", 'ignorecase', true )
    
    opt.fminopt = oc.fminopt;
    opt.gaopt = oc.gaopt;
    opt.gaopt.UseParallel = true;

    pp = gcp('nocreate');
    numcores = feature('numcores');
    if isempty(pp)
        parpool( min( [ oc.numWorkers, numcores ] ) );
    else
        if pp.NumWorkers < oc.numWorkers
            delete( pp );
            parpool( min( [ oc.numWorkers, numcores] ) );
        end
    end

    opt.popt = opt.generateProxy( opt, oc );
    [ oc, pulse, opt.popt ] = processAdjointFunctions( oc, pulse, opt.popt );
    opt.popt.generateWaveforms = opt.generateProxyWaveform;

end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ opt, oc ] = pSc0Process( opt, oc, ensureFeasibleStart )
if isfield( oc, "pSc0" ) || isfield( opt, "pSc0" )
    if isfield( oc, "pSc0" ) && ~isfield( opt, "pSc0" )
        opt.pSc0 = oc.pSc0;
    elseif isfield( opt, "pSc0" ) && isfield( oc, "pSc0" )
        oc.pSc0 = opt.pSc0;
    elseif isfield( opt, "pSc0" ) && ~isfield( oc, "pSc0" )
        oc.pSc0 = opt.pSc0;
    end

    opt.p0 = opt.scVec .*  opt.pSc0;

    if ensureFeasibleStart
        st = struct;
        if isfield( oc, 'ensureFeasibleStartRelLineSrchBnd' )
            st.relLineSrchBnd = oc.ensureFeasibleStartRelLineSrchBnd;
        end
        if isfield( oc, 'ensureFeasibleStartScaleAtEnd' )
            st.scaleAtEnd = oc.ensureFeasibleStartScaleAtEnd;
        end
        
        [ opt.pSc0, ~ ] = ensureFeasibleStartFMIN( opt.pSc0, opt, st );
        
        opt.p0 = opt.scVec .*  opt.pSc0;
    end

else
    % warning("Initializing pSc0 with zero vector of size:\t[%i,1]",...
    %     opt.numVars);
    % opt.pSc0 = zeros( opt.numVars, 1 );

    pScMax = 1e-4;
    warning("Initializing pSc0 with randn (mean: %g) vector of size:\t[%i,1]",...
        pScMax, opt.numVars);

    opt.pSc0 = pScMax * randn( opt.numVars, 1 );
    if ensureFeasibleStart
        [ opt.pSc0, ~ ] = ensureFeasibleStartFMIN( opt.pSc0, opt );
        opt.p0 = opt.scVec .*  opt.pSc0;
    end
end
end
% ----------------------------------------------------------------------- %