function [ RFOpt, resid, ompSt ] = solveOMPMatricesSpokes( AE, bE, RF0p, spokes, opt, curr_numSpokes, fminopt )

curr_bRF_idx = ( spokes.numSpokes * opt.numXYCoils - curr_numSpokes * opt.numXYCoils + 1 ):( spokes.numSpokes * opt.numXYCoils );

% change indices for the smaller kTP problem
breal_idx = transpose( 1 : (curr_numSpokes * opt.numXYCoils) ); % opt.breal_idx( curr_bRF_idx );
bimag_idx = transpose( ( curr_numSpokes * opt.numXYCoils + 1 ) : ( 2 * curr_numSpokes * opt.numXYCoils ) ); %opt.bimag_idx( curr_bRF_idx );
scVec = opt.scVec( [ opt.breal_idx( curr_bRF_idx ); opt.bimag_idx( curr_bRF_idx ) ] );
numVars = ( 2 * curr_numSpokes * opt.numXYCoils );

opt.breal_idx = breal_idx;
opt.bimag_idx = bimag_idx;
opt.numSpokes = curr_numSpokes;
opt.scVec = scVec;
opt.numVars = numVars;

% Update possible nonlinear constraints
if ~isempty( opt.nlconEqFuncs ) || ~isempty( opt.nlconIneqFuncs )
    [ nonlconIneqFuncsSpokes, nonlconEqFuncsSpokes ] = nlconAdjustSpokes( opt );
    opt.nonlcon = @( pSc ) nonlconSpokesCombine( pSc( : ), opt, spokes, nonlconIneqFuncsSpokes, nonlconEqFuncsSpokes );
else
    opt.nonlcon = [];
end

if nargin < 6
    fminopt = optimoptions( "fmincon" );
    fminopt.MaxFunctionEvaluations = inf;
    fminopt.Display = "off";
    fminopt.SpecifyConstraintGradient = true;
    fminopt.SpecifyObjectiveGradient = true;
    fminopt.StepTolerance = 1e-5;
    fminopt.FunctionTolerance = 1e-6;
    fminopt.ConstraintTolerance = 1e-14;
    fminopt.OptimalityTolerance = 1e-7;

    fminopt.MaxIterations = 500;

    fminopt.Algorithm = "sqp";

    % fminopt.Algorithm = "active-set";
    % fminopt.RelLineSrchBnd = 1e-1;
    % fminopt.RelLineSrchBndDuration = inf;
    % fminopt.TolConSQP = 1e-10;
    
end

Aopt = [];
bopt = [];
Aeqopt = [];
beqopt = [];
lbopt = -1 * ones( 2*opt.numSpokes*opt.numXYCoils, 1 );
ubopt = 1 * ones( 2*opt.numSpokes*opt.numXYCoils, 1 );
nonlconopt = opt.nonlcon;

costFnAnon = @( RFunsc ) residualCostFunction( RFunsc, AE, bE, opt );

RF0 = [...
    zeros( opt.numXYCoils, 1 ); real( RF0p );...
    zeros( opt.numXYCoils, 1 ); imag( RF0p ) ];
RFsc0 = RF0 ./ scVec;

ompSt = struct;
[ RFstackunscOpt, ompSt.fval, ompSt.exitflag, ompSt.output ] = fmincon( costFnAnon, RFsc0,...
    Aopt, bopt, Aeqopt, beqopt, lbopt, ubopt, nonlconopt, fminopt );

RFstackOpt = scVec .* RFstackunscOpt;

RFOpt = complex( RFstackOpt( breal_idx ), RFstackOpt( bimag_idx ) );

resid = AE * RFOpt - bE;
end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ cost, gradcost ] = residualCostFunction( RFunsc, AE, bE, opt )

breal_RF_idx = transpose( 1:(opt.numSpokes*opt.numXYCoils ) );
bimag_RF_idx = transpose( (opt.numSpokes*opt.numXYCoils + 1 ):(2*opt.numSpokes*opt.numXYCoils ) );
breal_sc = opt.scVec( opt.breal_idx );
bimag_sc = opt.scVec( opt.bimag_idx );

RFbrealunsc = RFunsc( breal_RF_idx );
RFbimagunsc = RFunsc( bimag_RF_idx );

RF = complex( breal_sc .* RFbrealunsc, bimag_sc .* RFbimagunsc );

Ax = AE * RF;
Axmb = Ax - bE;

norm2resid = norm( Axmb, 2 )^2;
norm2bspatial = norm( bE, 2 )^2;
cost = norm2resid / norm2bspatial;

if nargout > 1
    gradcost = zeros( 2*opt.numSpokes*opt.numXYCoils, 1 );
    gradcost_unsc = 2 * ( ctranspose( AE ) * ( AE*RF - bE ) );
    gradcost( breal_RF_idx ) = real( gradcost_unsc ) .* ( breal_sc / norm2bspatial );
    gradcost( bimag_RF_idx ) = imag( gradcost_unsc ) .* ( bimag_sc / norm2bspatial );
end

end
% ----------------------------------------------------------------------- %