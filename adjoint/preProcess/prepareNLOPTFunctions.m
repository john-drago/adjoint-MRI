function [ opt, nlopt ] = prepareNLOPTFunctions( opt, nlopt )

%% Handle bounds
nlopt.lower_bounds = opt.lb;
nlopt.upper_bounds = opt.ub;

%% Define objective function
min_objective = @( pSc ) objectiveFunction( pSc, opt );
nlopt.min_objective = min_objective;

if nlopt.algorithm == NLOPT_LD_SLSQP

    %% Define inequality constraints
    if ( ~isempty( opt.A ) ) || ( ~isempty( opt.nlconIneqFuncs ) )

        mfc = @( pSc ) constraintSLSQPFn( pSc, opt.A, opt.b, opt.nlconIneqFuncs, opt );
        nlopt.mfc = mfc;

        if ( ~isempty( opt.A ) )
            numAConstraints = size( opt.A, 1 );
        else
            numAConstraints = 0;
        end
        if ( ~isempty( opt.nlconIneqFuncs ) )
            numnlIneqConstraints = opt.nlconIneqCumSum( end );
        else
            numnlIneqConstraints = 0;
        end
        numIneqConstraints = numAConstraints + numnlIneqConstraints;

        nlopt.mfc_tol = nlopt.ConstraintTolerance * ones( numIneqConstraints, 1 );
        nlopt.mfc_count = numIneqConstraints;

    end

    %% Define equality constraints
    if ( ~isempty( opt.Aeq ) ) || ( ~isempty( opt.nlconEqFuncs ) )
        mh = @( pSc ) constraintSLSQPFn( pSc, opt.Aeq, opt.beq, opt.nlconEqFuncs, opt );

        nlopt.mh = mh;

        if ( ~isempty( opt.Aeq ) )
            numAeqConstraints = size( opt.Aeq, 1 );
        else
            numAeqConstraints = 0;
        end
        if ( ~isempty( opt.nlconEqFuncs ) )
            numnlEqConstraints = opt.nlconEqCumSum( end );
        else
            numnlEqConstraints = 0;
        end
        numEqConstraints = numAeqConstraints + numnlEqConstraints;

        nlopt.mh_tol = nlopt.ConstraintTolerance * ones( numEqConstraints, 1 );
        nlopt.mh_count = numEqConstraints;

    end

elseif nlopt.algorithm == NLOPT_LD_MMA

    %% Define inequality constraints
    if ( ~isempty( opt.A ) ) || ( ~isempty( opt.Aeq ) ) ...
            || ( ~isempty( opt.nlconIneqFuncs ) ) || ( ~isempty( opt.nlconEqFuncs ) )

        mfc = @( pSc ) constraintMMAFn( pSc, opt.A, opt.b, opt.Aeq, opt.beq, opt.nlconIneqFuncs, opt.nlconEqFuncs, opt );
        nlopt.mfc = mfc;

        if ( ~isempty( opt.A ) )
            numAConstraints = size( opt.A, 1 );
        else
            numAConstraints = 0;
        end
        if ( ~isempty( opt.Aeq ) )
            numAeqConstraints = 2 * size( opt.Aeq, 1 );
        else
            numAeqConstraints = 0;
        end
        if ( ~isempty( opt.nlconIneqFuncs ) )
            numnlIneqConstraints = opt.nlconIneqCumSum( end );
        else
            numnlIneqConstraints = 0;
        end
        if ( ~isempty( opt.nlconEqFuncs ) )
            numnlEqConstraints = 2 * opt.nlconEqCumSum( end );
        else
            numnlEqConstraints = 0;
        end

        numIneqConstraints = numAConstraints + numAeqConstraints +...
            numnlIneqConstraints + numnlEqConstraints;

        nlopt.mfc_tol = nlopt.ConstraintTolerance * ones( numIneqConstraints, 1 );
        nlopt.mfc_count = numIneqConstraints;

    end

else
    error( "Can't process this algorithm." )
end
end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ J, dJdp ] = objectiveFunction( pSc, opt )

costFn = opt.runCostFunction;

if nargout > 1
    [ J, gradpSc_J ] = costFn( pSc( : ), opt );
    
    dJdp = transpose( gradpSc_J );
    % dJdp = gradpSc_J;

    gradTol = 1e-10;
    dJdp( abs( dJdp ) < gradTol ) = 0;
    
else
    J = costFn( pSc( : ), opt );
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ c, gradc ] = constraintMMAFn( pSc, A, b, Aeq, beq, nlconIneqFuncs, nlconEqFuncs, opt )

pSc = pSc( : );

if ~isempty( A )
    cA = A * pSc - b;
    if nargout > 1
        gradcA = A;
    end
else
    cA = [];
    gradcA = [];
end

if ~isempty( Aeq )
    cAeq = [ ( Aeq * pSc - beq ) ; -( Aeq * pSc - beq ) ];
    if nargout > 1
        gradcAeq = [ Aeq; -Aeq ];
    end
else
    cAeq = [];
    gradcAeq = [];
end

cnlIneq = [];
gradcnlIneq = [];

if ~isempty( nlconIneqFuncs )

    for cc = 1:length( nlconIneqFuncs )
        nlconIneqFunc = nlconIneqFuncs{ cc };
        if nargout > 1
            [ cnlIneqi, gradcnlIneqi ] = nlconIneqFunc( pSc, opt );

            gradcnlIneq = [ gradcnlIneq; transpose( gradcnlIneqi ) ]; %#ok
        else
            cnlIneqi = nlconIneqFunc( pSc, opt );
        end
        cnlIneq = [ cnlIneq; cnlIneqi ]; %#ok
    end
end

cnlEq = [];
gradcnlEq = [];

if ~isempty( nlconEqFuncs )

    for cc = 1:length( nlconEqFuncs )
        nlconEqFunc = nlconEqFuncs{ cc };
        if nargout > 1
            [ cnlEqi, gradcnlEqi ] = nlconEqFunc( pSc, opt );

            gradcnlEq = [ gradcnlEq; transpose( gradcnlEqi ); -transpose( gradcnlEqi ) ]; %#ok
        else
            cnlEqi = nlconEqFunc( pSc, opt );
        end
        cnlEq = [ cnlEq; cnlEqi; -cnlEqi ]; %#ok
    end
end

c = [ cA; cAeq; cnlIneq; cnlEq ];

if nargout > 1
    gradc = [ gradcA; gradcAeq; gradcnlIneq; gradcnlEq ];

    gradTol = 1e-10;
    gradc( abs( gradc ) < gradTol ) = 0;
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ c, gradc ] = constraintSLSQPFn( pSc, A, b, nlconFuncs, opt )

pSc = pSc( : );

if ~isempty( A )
    cA = A * pSc - b;
    if nargout > 1
        gradcA = A;
    end
else
    cA = [];
    gradcA = [];
end

cnl = [];
gradcnl = [];
if ~isempty( nlconFuncs )

    for cc = 1:length( nlconFuncs )
        nlconFunc = nlconFuncs{ cc };
        if nargout > 1
            [ cnli, gradcnli ] = nlconFunc( pSc, opt );

            gradcnl = [ gradcnl; transpose( gradcnli ) ]; %#ok
        else
            cnli = nlconFunc( pSc, opt );
        end
        cnl = [ cnl; cnli ]; %#ok
    end
end

c = [ cA; cnl ];
if nargout > 1
    gradc = [ gradcA; gradcnl ];

    gradTol = 1e-10;
    gradc( abs( gradc ) < gradTol ) = 0;
end

end
% ----------------------------------------------------------------------- %