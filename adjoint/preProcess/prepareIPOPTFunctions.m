function [ opt, ipopt ] = prepareIPOPTFunctions( opt, ipopt )

%% Determine whether to output to file
% Determine if there is a savedir already
% timestr = strcat('__', char(datetime('now', 'format','yyMMddHHmmss')) );
if isfield( opt, 'saveDir' )
    saveDirTemp = opt.saveDir;
end
if opt.saveResult
    if ~isfield( ipopt.options.ipopt, 'output_file' )
        ipopt.options.ipopt.output_file = char( fullfile( saveDirTemp, "OUTPUT.txt" ) );
    end
    ipopt.options.ipopt.file_print_level = 5;
    if ~isfolder( saveDirTemp )
        mkdir( saveDirTemp )
    end
end

ipopt.funcs = struct;

%% Make sparse matrices
opt.A = sparse( opt.A );
opt.Aeq = sparse( opt.Aeq );

%% Define constraints function
% Define constraints constraints
if ( ~isempty( opt.A ) ) || ( ~isempty( opt.Aeq ) ) ...
        || ( ~isempty( opt.nlconIneqFuncs ) ) || ( ~isempty( opt.nlconEqFuncs ) )

    constraints = @( pSc ) constraintFn( pSc, opt.A, opt.b, opt.Aeq, opt.beq, opt.nlconIneqFuncs, opt.nlconEqFuncs, opt );
    ipopt.funcs.constraints = constraints;

    if ( ~isempty( opt.A ) )
        numAConstraints = size( opt.A, 1 );
    else
        numAConstraints = 0;
    end
    if ( ~isempty( opt.Aeq ) )
        numAeqConstraints = size( opt.Aeq, 1 );
    else
        numAeqConstraints = 0;
    end
    if ( ~isempty( opt.nlconIneqFuncs ) )
        numnlIneqConstraints = opt.nlconIneqCumSum( end );
    else
        numnlIneqConstraints = 0;
    end
    if ( ~isempty( opt.nlconEqFuncs ) )
        numnlEqConstraints = opt.nlconEqCumSum( end );
    else
        numnlEqConstraints = 0;
    end

    numConstraints = numAConstraints + numAeqConstraints +...
        numnlIneqConstraints + numnlEqConstraints;
    
    ipopt.options.cl = -inf * ones( numConstraints, 1 );
    ipopt.options.cu = inf * ones( numConstraints, 1 );
    
    % Define lower bounds on constraints
    ipopt.options.cl( (numAConstraints+1):...
        (numAConstraints+numAeqConstraints), 1 ) = 0;
    ipopt.options.cl( (numAConstraints+numAeqConstraints+numnlIneqConstraints+1):...
        (numAConstraints+numAeqConstraints+numnlIneqConstraints+numnlEqConstraints), 1 ) = 0;
    
    % Define upper bounds on constraints
    ipopt.options.cu( 1:numAConstraints, 1 ) = 0;
    ipopt.options.cu( (numAConstraints+1):...
        (numAConstraints+numAeqConstraints), 1 ) = 0;
    ipopt.options.cu( (numAConstraints+numAeqConstraints+1):...
        (numAConstraints+numAeqConstraints+numnlIneqConstraints), 1 ) = 0;
    ipopt.options.cu( (numAConstraints+numAeqConstraints+numnlIneqConstraints+1):...
        (numAConstraints+numAeqConstraints+numnlIneqConstraints+numnlEqConstraints), 1 ) = 0;

end

%% Add constraints opt information
[ ipopt.funcs.constraintNames, ipopt.funcs.constraintAmts, ipopt.funcs.constraintCumSum ] = ...
    getIPOPTConstraintLocations( opt );

%% Define jacobian of constraints function
ipopt.funcs.jacobian = @( pSc ) constraintJacobianFn( pSc, opt.A, opt.Aeq, opt.nlconIneqFuncs, opt.nlconEqFuncs, opt );

%% Define jacobian structure
jsCenterVal = 0.5;
jsDevVal = 0.25;
jsVals = jsCenterVal + ( 2 * jsDevVal * rand( opt.numVars, 1 ) - jsDevVal ) - eps( 1e5 );
jsVals = jsVals .* ( 2 * double( randi( [0, 1], size( jsVals ) ) ) - 1 );

jacobianStructure = sparse( abs( ipopt.funcs.jacobian( jsVals ) ) > 0 );

if ~( numAConstraints == 0 ) || ~( numAeqConstraints == 0 )
    jacobianStructure( numAConstraints+numAeqConstraints, : ) = true;
end

%% Find variable for RF indices
RFSARvarIdxsLog = matches( opt.varNames,...
    [ "breal"; "bimag"; "b"; "breal-ORSP"; "bimag-ORSP"; "breal-MPSP"; "bimag-MPSP"; "dt" ],...
    "ignorecase", true );
RFSARvarIdxs = cell2mat( opt.varIdxs( RFSARvarIdxsLog, 1 ) );

% Find constraints for RF, SAR, and power
RFSARconstraintID = find( ...
    contains( string( ipopt.funcs.constraintNames ),...
    [ "RFPower"; "SAR"; "RFMax"; "MaxRF"; "RFSlew"; "RFAccel"; "RFBandwidth" ], 'IgnoreCase', true ) );

if any( RFSARconstraintID == 1 )
    RFSARconstIdx = 1:ipopt.funcs.constraintAmts( RFSARconstraintID( 1 ) );
else
    RFSARconstIdx = [];
end

RFSARconstraintIDnot1 = RFSARconstraintID( RFSARconstraintID ~= 1 );

for ii = 1 : length( RFSARconstraintIDnot1 )
    RFSARconstIdxii = (ipopt.funcs.constraintCumSum( RFSARconstraintIDnot1( ii )-1 )+1) : ( ipopt.funcs.constraintCumSum( RFSARconstraintIDnot1( ii )-1 ) +...
        ipopt.funcs.constraintAmts( RFSARconstraintIDnot1( ii ) ) );
    RFSARconstIdx = [ RFSARconstIdx; RFSARconstIdxii( : ) ]; %#ok
end

if ~isempty( RFSARconstraintID )
    RFSARconstraintIDMat = repmat( RFSARconstIdx, [ 1, length( RFSARvarIdxs )] );
    RFSARvarIdxsMat = repmat( transpose( RFSARvarIdxs ), [ length( RFSARconstIdx ), 1] );

    RFSARconstraintJSIdxs = sub2ind( size( jacobianStructure ), RFSARconstraintIDMat, RFSARvarIdxsMat );
    jacobianStructure( RFSARconstraintJSIdxs ) = true( size( RFSARconstraintJSIdxs ) );
end

%% Find variable for grad indices
gradVarIdxsLog = matches( opt.varNames,...
    [ "grad";...
    "Gx"; "Gy"; "Gz";...
    "Gx-ORSP"; "Gy-ORSP"; "Gz-ORSP";...
    "Gx-Blip"; "Gy-Blip"; "Gz-Blip";...
    "Gxreal-MPSP"; "Gyreal-MPSP"; "Gzreal-MPSP"; "Gximag-MPSP"; "Gyimag-MPSP"; "Gzimag-MPSP";...
    "kmax"; "a"; "b"; "u"; "v" ],...
    "ignorecase", true );

gradVarIdxs = cell2mat( opt.varIdxs( gradVarIdxsLog, 1 ) );

% Find constraints for RF, SAR, and power
gradConstraintID = find( ...
    contains( string( ipopt.funcs.constraintNames ),...
    [ "GradMax"; "MaxGrad"; "GradSlew"; "GradAccel";  ], 'IgnoreCase', true ) );

if any( gradConstraintID == 1 )
    gradConstIdx = 1:ipopt.funcs.constraintAmts( gradConstraintID( 1 ) );
else
    gradConstIdx = [];
end

gradConstraintIDnot1 = gradConstraintID( gradConstraintID ~= 1 );

for ii = 1 : length( gradConstraintIDnot1 )
    gradConstIdxii = (ipopt.funcs.constraintCumSum( gradConstraintIDnot1( ii )-1 )+1) : ( ipopt.funcs.constraintCumSum( gradConstraintIDnot1( ii )-1 ) +...
        ipopt.funcs.constraintAmts( gradConstraintIDnot1( ii ) ) );
    gradConstIdx = [ gradConstIdx; gradConstIdxii( : ) ]; %#ok
end

if ~isempty( gradConstraintID )
    gradConstraintIDMat = repmat( gradConstIdx, [ 1, length( gradVarIdxs )] );
    gradVarIdxsMat = repmat( transpose( gradVarIdxs ), [ length( gradConstIdx ), 1] );

    gradConstraintJSIdxs = sub2ind( size( jacobianStructure ), gradConstraintIDMat, gradVarIdxsMat );
    jacobianStructure( gradConstraintJSIdxs ) = true( size( gradConstraintJSIdxs ) );
end

%% Find variable for shim indices
shimVarIdxsLog = matches( opt.varNames,...
    [ "shim";...
    "shim-Blip";...
    "shimreal-MPSP"; "shimimag-MPSP";...
    "shimmag-MPSP"; "shimph-MPSP"; "shimblipscale" ],...
    "ignorecase", true );
shimVarIdxs = cell2mat( opt.varIdxs( shimVarIdxsLog, 1 ) );

% Find constraints for RF, SAR, and power
shimConstraintID = find( ...
    contains( string( ipopt.funcs.constraintNames ),...
    [ "ShimMax"; "MaxShim"; "ShimSlew"; "ShimAccel"; "ShimTotal" ], 'IgnoreCase', true ) );

if any( shimConstraintID == 1 )
    shimConstIdx = 1:ipopt.funcs.constraintAmts( shimConstraintID( 1 ) );
else
    shimConstIdx = [];
end

shimConstraintIDnot1 = shimConstraintID( shimConstraintID ~= 1 );

for ii = 1 : length( shimConstraintIDnot1 )
    shimConstIdxii = (ipopt.funcs.constraintCumSum( shimConstraintIDnot1( ii )-1 )+1) : ( ipopt.funcs.constraintCumSum( shimConstraintIDnot1( ii )-1 ) +...
        ipopt.funcs.constraintAmts( shimConstraintIDnot1( ii ) ) );
    shimConstIdx = [ shimConstIdx; shimConstIdxii( : ) ]; %#ok
end

if ~isempty( shimConstraintID )
    shimConstraintIDMat = repmat( shimConstIdx, [ 1, length( shimVarIdxs )] );
    shimVarIdxsMat = repmat( transpose( shimVarIdxs ), [ length( shimConstIdx ), 1] );

    shimConstraintJSIdxs = sub2ind( size( jacobianStructure ), shimConstraintIDMat, shimVarIdxsMat );
    jacobianStructure( shimConstraintJSIdxs ) = true( size( shimConstraintJSIdxs ) );
end

%% Assign jacobian structure
ipopt.funcs.jacobianstructure = @() jacobianStructureFn( jacobianStructure );

%% Define lower and upper bounds on constraints
ipopt.options.lb = opt.lb;
ipopt.options.ub = opt.ub;

%% Define hessian approximation
ipopt.options.ipopt.hessian_approximation = 'limited-memory';
ipopt.options.ipopt.limited_memory_update_type = 'bfgs';

%% Define objective function
objectiveFn = @( pSc ) objectiveFunction( pSc, opt );
ipopt.funcs.objective = @( pSc ) ipoptObjectiveBestOptWrapper( pSc, objectiveFn, ipopt.funcs, ipopt.options );

%% Define gradient
ipopt.funcs.gradient = @( pSc ) gradientFunction( pSc, opt );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ constraintNames, constraintAmts, constraintCumSum ] = ...
    getIPOPTConstraintLocations( opt )

constraintNames = strings( 0, 0 );
constraintAmts = [];
constraintIdx = 0;

% Ab inequality
if ~isempty( opt.AbConstraintNames )
    num_AbConstraints = size( opt.AbConstraintAmts, 1 );

    for nn = 1:num_AbConstraints
        constraintIdx = constraintIdx + 1;
        constraintAmts( constraintIdx, 1 ) = uint32( opt.AbConstraintAmts( nn, 1 ) ); %#ok
        constraintNames( constraintIdx, 1 ) = sprintf( "Ab-%s", opt.AbConstraintNames( nn, 1 ) );
    end
end

% Ab equality
if ~isempty( opt.AbeqConstraintNames )
    num_AbeqConstraints = size( opt.AbeqConstraintAmts, 1 );

    for nn = 1:num_AbeqConstraints
        constraintIdx = constraintIdx + 1;
        constraintAmts( constraintIdx, 1 ) = uint32( opt.AbeqConstraintAmts( nn, 1 ) ); %#ok
        constraintNames( constraintIdx, 1 ) = sprintf( "Abeq-%s", opt.AbeqConstraintNames( nn, 1 ) );
    end
end

% nlcon inequality
if ~isempty( opt.nlconIneqFuncNames )
    num_nlconIneqConstraints = size( opt.nlconIneqAmts, 1 );

    for nn = 1:num_nlconIneqConstraints
        constraintIdx = constraintIdx + 1;
        constraintAmts( constraintIdx, 1 ) = uint32( opt.nlconIneqAmts( nn, 1 ) ); %#ok
        constraintNames( constraintIdx, 1 ) = sprintf( "nlconIneq-%s", opt.nlconIneqFuncNames( nn, 1 ) );
    end
end

% nlcon equality
if ~isempty( opt.nlconEqFuncNames )
    num_nlconEqConstraints = size( opt.nlconEqAmts, 1 );

    for nn = 1:num_nlconEqConstraints
        constraintIdx = constraintIdx + 1;
        constraintAmts( constraintIdx, 1 ) = uint32( opt.nlconEqAmts( nn, 1 ) ); %#ok
        constraintNames( constraintIdx, 1 ) = sprintf( "nlconEq-%s", opt.nlconEqFuncNames( nn, 1 ) );
    end
end

constraintCumSum = cumsum( constraintAmts );

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function jacobianStructure = jacobianStructureFn( jacobianStructure )
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ gradc ] = constraintJacobianFn( pSc, A, Aeq, nlconIneqFuncs, nlconEqFuncs, opt )

pSc = pSc( : );

if ~isempty( A )
    gradcA = A;
else
    gradcA = [];
end

if ~isempty( Aeq )
    gradcAeq = Aeq;
else
    gradcAeq = [];
end

gradcnlIneq = [];

if ~isempty( nlconIneqFuncs )

    for cc = 1:length( nlconIneqFuncs )
        nlconIneqFunc = nlconIneqFuncs{ cc };
        [ ~, gradcnlIneqi ] = nlconIneqFunc( pSc, opt );
        gradcnlIneq = [ gradcnlIneq; sparse( transpose( gradcnlIneqi ) ) ]; %#ok
    end

end

gradcnlEq = [];

if ~isempty( nlconEqFuncs )

    for cc = 1:length( nlconEqFuncs )
        nlconEqFunc = nlconEqFuncs{ cc };
        [ ~, gradcnlEqi ] = nlconEqFunc( pSc, opt );
        gradcnlEq = [ gradcnlEq; sparse( transpose( gradcnlEqi ) ) ]; %#ok
    end

end

gradc = [ gradcA; gradcAeq; gradcnlIneq; gradcnlEq ];

gradTol = 1e-10;
gradc( abs( gradc ) < gradTol ) = 0;

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ c ] = constraintFn( pSc, A, b, Aeq, beq, nlconIneqFuncs, nlconEqFuncs, opt )

pSc = pSc( : );

if ~isempty( A )
    cA = A * pSc - b;
else
    cA = [];
end

if ~isempty( Aeq )
    cAeq = ( Aeq * pSc - beq );
else
    cAeq = [];
end

cnlIneq = [];

if ~isempty( nlconIneqFuncs )
    for cc = 1:length( nlconIneqFuncs )
        nlconIneqFunc = nlconIneqFuncs{ cc };
        cnlIneqi = nlconIneqFunc( pSc, opt );
        cnlIneq = [ cnlIneq; cnlIneqi ]; %#ok
    end
end

cnlEq = [];

if ~isempty( nlconEqFuncs )

    for cc = 1:length( nlconEqFuncs )
        nlconEqFunc = nlconEqFuncs{ cc };
        cnlEqi = nlconEqFunc( pSc, opt );
        cnlEq = [ cnlEq; cnlEqi ]; %#ok
    end
end

c = [ cA; cAeq; cnlIneq; cnlEq ];

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ gradJ_pSc ] = gradientFunction( pSc, opt )

% [ ~, gradJ_pSc ] = objectiveFn( pSc( : ), opt );

costFn = opt.runCostFunction;
[ ~, gradJ_pSc ] = costFn( pSc( : ), opt );

gradTol = 1e-10;
gradJ_pSc( abs( gradJ_pSc ) < gradTol ) = 0;

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ J, gradJ_pSc ] = objectiveFunction( pSc, opt )

% if nargout > 1
%     [ J, gradJ_pSc ] = opt.runCostFunction( pSc( : ), opt );
% else
%     J = opt.runCostFunction( pSc( : ), opt );
% end

costFn = opt.runCostFunction;
if nargout > 1
    [ J, gradJ_pSc ] = costFn( pSc( : ), opt );

    gradTol = 1e-10;
    gradJ_pSc( abs( gradJ_pSc ) < gradTol ) = 0;
else
    J = costFn( pSc( : ), opt );
end

end
% ----------------------------------------------------------------------- %