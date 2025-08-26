function [ RFopt, Gopt, Kopt, Aspatial, bspatial ] = optimizekTPointsSTA(...
    RF0, K0, RFMax, gradSlewRate, opt, pulse, maxIterations )

if nargin < 7
    maxIterations = 500;
end

%% Initialization
num_kTP = pulse.num_kTP;

%% Specify Gradient Optimization Options
fminopt = optimoptions( "fmincon" );
fminopt.MaxFunctionEvaluations = inf;
fminopt.FunctionTolerance = 1e-8;
fminopt.StepTolerance = 1e-8;
fminopt.ConstraintTolerance = 1e-9;
fminopt.SpecifyObjectiveGradient = true;
fminopt.SpecifyConstraintGradient = true;
fminopt.Display = "iter";

fminopt.MaxIterations = maxIterations;

% fminopt.Algorithm = "sqp";

fminopt.Algorithm = "active-set";
fminopt.RelLineSrchBnd = 1e-2;
fminopt.RelLineSrchBndDuration = inf;
fminopt.TolConSQP = 1e-10;

%% Get number of variables
numVars = ( num_kTP - 1 ) * 3 + 2 * num_kTP * opt.numXYCoils;

%% Set indices
staSt = struct;
k_idx = uint32( transpose( 1: ( ( num_kTP - 1 ) * 3 ) ) );
breal_idx = uint32( transpose( ( ( num_kTP - 1 ) * 3 + 1 ) : ( ( num_kTP - 1 ) * 3 + num_kTP * opt.numXYCoils ) ) );
bimag_idx = uint32( transpose( ( ( ( num_kTP - 1 ) * 3 + num_kTP * opt.numXYCoils + 1 ) ) : ( ( num_kTP - 1 ) * 3 + 2 * num_kTP * opt.numXYCoils ) ) );

staSt.k_idx = k_idx;
staSt.breal_idx = breal_idx;
staSt.bimag_idx = bimag_idx;

%% Specify inequalities
lbopt = -inf * ones( numVars, 1 );
ubopt = inf * ones( numVars, 1 );
scVec = zeros( numVars, 1 );

lbopt( breal_idx ) = - ones( length( breal_idx ), 1 );
ubopt( breal_idx ) = ones( length( breal_idx ), 1 );
scVec( breal_idx ) = RFMax / sqrt( 2 ) * ones( length( breal_idx ), 1 );

lbopt( bimag_idx ) = - ones( length( bimag_idx ), 1 );
ubopt( bimag_idx ) = ones( length( bimag_idx ), 1 );
scVec( bimag_idx ) = RFMax / sqrt( 2 ) * ones( length( bimag_idx ), 1 );

kScale = 50;
scVec( k_idx ) = kScale;

staSt.scVec = scVec;

%% Try to incorporate nonlinear constraints for the pulse design
opt.breal_idx = breal_idx;
opt.bimag_idx = bimag_idx;
opt.scVec = scVec;
opt.numVars = numVars;

% Update possible nonlinear constraints  
if ~isempty( opt.nlconEqFuncs ) || ~isempty( opt.nlconIneqFuncs )
    opt.nonlcon = @( pSc ) nonlconCombine( pSc( : ), opt );
else
    opt.nonlcon = [];
end

%% Specify Constraints
if num_kTP > 1
    Apos_k =...
        diag( kScale * ones( 3 * ( num_kTP - 1 ), 1 ), 0 ) +...
        diag( ( kScale * -1 ) * ones( 3 * ( num_kTP - 2 ), 1 ), 3 );
    bpos_k = ( ( gradSlewRate * pulse.blipLength^2 * opt.gyro ) / ( 8 * pi ) ) * ones( 3 * ( num_kTP - 1 ), 1 );

    Apos_k = [ Apos_k, zeros( size( Apos_k, 1 ), 2 * num_kTP * opt.numXYCoils ) ];

    A_k = [ Apos_k; -Apos_k ];
    b_k = [ bpos_k; bpos_k ];
else
    A_k = [];
    b_k = [];
end

Aopt = A_k;
bopt = b_k;
Aeqopt = [];
beqopt = [];
nonlconopt = opt.nonlcon;

%% Initial guess
K0vec = K0( :, 1:(end-1) );
K0vec = K0vec( : );

RF0vec = RF0( : );

sta0 = zeros( numVars, 1 );
sta0( k_idx ) = K0vec;
sta0( breal_idx ) = real( RF0vec );
sta0( bimag_idx ) = imag( RF0vec );
staSc0 = sta0 ./ scVec;

%% Specify Objective Function
MLSObjFunAnon = @( sta ) MLSObjFun( sta, staSt, opt, pulse );

%% Run optimization
[ staScopt, ~, ~, ~ ] = fmincon( MLSObjFunAnon, staSc0, Aopt, bopt, Aeqopt, beqopt, lbopt, ubopt, nonlconopt, fminopt );
% [ staScopt, cost, exitflag, output ] = fmincon( MLSObjFunAnon, staSc0, Aopt, bopt, Aeqopt, beqopt, lbopt, ubopt, nonlconopt, fminopt );
staopt = staSt.scVec .* staScopt;

%% Rearrange k-space vector
Kopt = [ reshape( staopt( staSt.k_idx ), [ 3, (num_kTP - 1) ] ), zeros( 3, 1 ) ];

%% Get Gradients
Gopt = calculateGradientBlipskTPSTA( Kopt, opt, pulse );

%% Get RF waveforms
RFoptvec = complex( staopt( staSt.breal_idx ), staopt( staSt.bimag_idx ) );
RFopt = reshape( RFoptvec, [ opt.numXYCoils, num_kTP ] );

%% Remake optimization matrices
if nargout > 3
    [ Aspatial, bspatial ] = generatekTPSTAMatrices( Kopt, pulse, opt );
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ cost, gradcost ] = MLSObjFun( staSc, staSt, opt, pulse )

num_kTP = pulse.num_kTP;

Kvec = staSt.scVec( staSt.k_idx ) .* staSc( staSt.k_idx );
K = [ reshape( Kvec, [ 3, ( pulse.num_kTP - 1 ) ] ), zeros( 3, 1 ) ];

%% Prepare Matrices
[ Aspatial, bspatial ] = generatekTPSTAMatrices( K, pulse, opt );

RFrealvec = staSt.scVec( staSt.breal_idx ) .* staSc( staSt.breal_idx );
RFimagvec = staSt.scVec( staSt.bimag_idx ) .* staSc( staSt.bimag_idx );
RFcompvec = complex( RFrealvec, RFimagvec );
% RF = reshape( RFcompvec, [ opt.numXYCoils, pulse.num_kTP ] );

Ax = Aspatial * RFcompvec;
absAx = abs( Ax );
magb = abs( bspatial );
absAxmb = absAx - magb;

norm2resid = norm( absAxmb, 2 )^2;
norm2bspatial = norm( magb, 2 )^2;
cost = norm2resid / norm2bspatial;

if nargout > 1

    % pre-calculate values
    dfdabsAx = 2 * transpose( absAxmb .* ( conj( Ax ) ./ absAx ) );
    
    % deal with RF 
    dfdRFan_comp = transpose( conj( dfdabsAx * Aspatial ) );
    dfdRFan_real = real( dfdRFan_comp );
    dfdRFan_imag = imag( dfdRFan_comp );

    % deal with grad
    Xadj = ( opt.pos( :, 1 ) * 2 * pi ) * 1j;
    Yadj = ( opt.pos( :, 2 ) * 2 * pi ) * 1j;
    Zadj = ( opt.pos( :, 3 ) * 2 * pi ) * 1j;

    dAxdkx = zeros( opt.numPos, (num_kTP-1) );
    dAxdky = zeros( opt.numPos, (num_kTP-1) );
    dAxdkz = zeros( opt.numPos, (num_kTP-1) );

    RFidx = reshape( 1:(opt.numXYCoils * num_kTP ) , [ opt.numXYCoils, num_kTP ] );

    % for nn = 1:(num_kTP-1)
    %     dAxdkx( :, nn ) = ( Xadj .* Aspatial( :, RFidx( :, nn ) ) ) * RFcompvec( RFidx( :, nn ), 1 );
    %     dAxdky( :, nn ) = ( Yadj .* Aspatial( :, RFidx( :, nn ) ) ) * RFcompvec( RFidx( :, nn ), 1 );
    %     dAxdkz( :, nn ) = ( Zadj .* Aspatial( :, RFidx( :, nn ) ) ) * RFcompvec( RFidx( :, nn ), 1 );
    % end

    for nn = 1:(num_kTP-1)
        dAxdkx( :, nn ) = Xadj .* ( Aspatial( :, RFidx( :, nn ) ) * RFcompvec( RFidx( :, nn ), 1 ) );
        dAxdky( :, nn ) = Yadj .* ( Aspatial( :, RFidx( :, nn ) ) * RFcompvec( RFidx( :, nn ), 1 ) );
        dAxdkz( :, nn ) = Zadj .* ( Aspatial( :, RFidx( :, nn ) ) * RFcompvec( RFidx( :, nn ), 1 ) );
    end
    
    dfdkanrshp = [...
        real( dfdabsAx * dAxdkx );...
        real( dfdabsAx * dAxdky );...
        real( dfdabsAx * dAxdkz );...
        ];

    dfdkan = dfdkanrshp( : );

    gradcost_unsc = zeros( size( staSc, 1 ), 1 );
    gradcost_unsc( staSt.breal_idx ) = dfdRFan_real .* staSt.scVec( staSt.breal_idx );
    gradcost_unsc( staSt.bimag_idx ) = dfdRFan_imag .* staSt.scVec( staSt.bimag_idx );
    gradcost_unsc( staSt.k_idx ) = dfdkan .* staSt.scVec( staSt.k_idx );
    gradcost = gradcost_unsc / norm2bspatial;

end

end
% ----------------------------------------------------------------------- %
