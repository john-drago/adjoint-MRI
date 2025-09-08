function [ RFopt, Gopt, Kopt, sopt, staSt, output ] = optimizeSPINSSTA(...
    initParms, RFMax, gradSlewRate, opt, pulse, initialFminMaxIterations, finalFminMaxIterations )

%% Initialize Function
if nargin < 6
    initialFminMaxIterations = 100;
end
if nargin < 7
    finalFminMaxIterations = 200;
end

opt.gradSlewRate_constr = gradSlewRate;

%% Create tvec
timeResDigits = 6; % round to microseconds
dt_tol = 1e-6;
pulseLength = round( pulse.length, timeResDigits); 
dt = opt.dt;

% Determine if time point spacing is a good interval
if abs( ( round( pulseLength/dt, 10 ) - pulseLength/dt ) / pulseLength ) < dt_tol
    tvec = (dt/2) : dt : pulseLength;
else
    dt = pulseLength / round( pulseLength/dt );
    tvec = (dt/2) : dt : pulseLength;
end

dtvec = dt * ones( size( tvec ) );
numTimePoints = length( tvec );

%% Handle SPINS specific constraints
% kmax
if isfield( pulse, "min_kmax" )
    min_kmax = pulse.min_kmax;
else
    min_kmax = 5; % rad/m
end
if isfield( pulse, "max_kmax" )
    max_kmax = pulse.max_kmax;
else
    max_kmax = 100; % rad/m
end

% u
if isfield( pulse, "min_u" )
    min_u = pulse.min_u;
else
    min_u = ( 2*pi ) / pulse.length; % rad/m
end
if isfield( pulse, "max_u" )
    max_u = pulse.max_u;
else
    max_u = ( 2*pi * 5 ) / pulse.length; % rad/m
end

% v
if isfield( pulse, "min_v" )
    min_v = pulse.min_v;
else
    min_v = ( 2*pi ) / pulse.length; % rad/m
end
if isfield( pulse, "max_v" )
    max_v = pulse.max_v;
else
    max_v = ( 2*pi * 5 ) / pulse.length; % rad/m
end

% a
if isfield( pulse, "min_a" )
    min_a = pulse.min_a;
else
    min_a = 0.5;
end
if isfield( pulse, "max_a" )
    max_a = pulse.max_a;
else
    max_a = 50;
end

% b
if isfield( pulse, "min_b" )
    min_b = pulse.min_b;
else
    min_b = 0.0;
end
if isfield( pulse, "max_b" )
    max_b = pulse.max_b;
else
    max_b = 1.0;
end

%% Set optimization indices for initial
varArray = {...
        "breal", numTimePoints * opt.numXYCoils;...
        "bimag", numTimePoints * opt.numXYCoils;...
        };

staSt = struct;

[ staSt, ~, ~, ~, ~, ~, ~, ~ ] =...
        processVarOrganization( staSt, varArray );

staSt.dtvec = dtvec;
staSt.tvec = tvec;
staSt.dt = dt;

staSt.numXYCoils = opt.numXYCoils;
staSt.numZCoils = opt.numZCoils;
staSt.gyro = opt.gyro;
staSt.Z0 = opt.Z0;
staSt.dutyCycle = opt.dutyCycle;
if isfield( opt, 'VOPs' )
    staSt.VOPs = opt.VOPs;
    staSt.numVOPs = opt.numVOPs;
end
if isfield( opt, 'QGlobal' )
    staSt.QGlobal = opt.QGlobal;
end

%% Determine indices of different periods
initSlewTime = pulse.initSlewTime;
finalSlewTime = pulse.finalSlewTime;

idxtol = dt/100;
init_slew_i = find( ( tvec ) >= ( 0 - idxtol ) , 1, 'first' );
init_slew_f = find( ( tvec ) <= ( initSlewTime + idxtol ), 1, 'last' );
spins_i = find( ( tvec ) >= ( initSlewTime - idxtol ), 1, 'first' );
spins_f = find( ( tvec ) <= ( pulseLength + idxtol ), 1, 'last' );
final_slew_i = find( ( tvec ) >= ( pulseLength - finalSlewTime - idxtol ), 1, 'first' );
final_slew_f = find( ( tvec ) <= ( pulseLength + idxtol ), 1, 'last' );

if init_slew_f >= spins_i
    spins_i = spins_i + 1;
end

% init_slew scaling
init_slew_sc = tvec( init_slew_i:init_slew_f ) / initSlewTime;
init_slew_num = init_slew_f - init_slew_i + 1;

% final_slew scaling
final_slew_sc = - ( tvec( final_slew_i:final_slew_f ) - pulseLength ) / finalSlewTime;
final_slew_num = final_slew_f - final_slew_i + 1;

% slew timing
spins_num = spins_f - spins_i + 1;
spins_int_i = spins_i;
spins_int_f = (final_slew_i - 1);
spins_int_idx = uint32( spins_int_i:spins_int_f );
spins_int_num = (final_slew_i - 1) - spins_i + 1;
Tspins = pulseLength - initSlewTime;

%% Assign to sta struct
staSt.initSlewTime = initSlewTime;
staSt.init_slew_i = init_slew_i;
staSt.init_slew_f = init_slew_f;
staSt.init_slew_sc = init_slew_sc;
staSt.init_slew_num = init_slew_num;

staSt.finalSlewTime = finalSlewTime;
staSt.final_slew_i = final_slew_i;
staSt.final_slew_f = final_slew_f;
staSt.final_slew_sc = final_slew_sc;
staSt.final_slew_num = final_slew_num;

staSt.spins_i = spins_i;
staSt.spins_f = spins_f;
staSt.spins_int_i = spins_int_i;
staSt.spins_int_f = spins_int_f;
staSt.spins_int_idx = spins_int_idx;
staSt.spins_int_num = spins_int_num;
staSt.spins_num = spins_num;

staSt.Tspins = Tspins;

staSt.numTimePoints = numTimePoints;
staSt.pulseLength = pulseLength;
staSt.tvec = tvec;
staSt.dtvec = dtvec;


%% Specify Gradient Optimization Options
fminopt = optimoptions( "fmincon" );
fminopt.MaxFunctionEvaluations = inf;
fminopt.StepTolerance = 1e-5;
fminopt.ConstraintTolerance = 1e-9;
fminopt.SpecifyObjectiveGradient = true;
fminopt.SpecifyConstraintGradient = true;

fminopt.Display = "iter";

% fminopt.Algorithm = "sqp";
fminopt.MaxIterations = initialFminMaxIterations;

fminopt.Algorithm = "active-set";
fminopt.RelLineSrchBnd = 1e-2;
fminopt.RelLineSrchBndDuration = inf;
fminopt.TolConSQP = 1e-10;

%% Specify inequalities
staSt.lb = -inf * ones( staSt.numVars, 1 );
staSt.ub = inf * ones( staSt.numVars, 1 );
staSt.scVec = zeros( staSt.numVars, 1 );

staSt.lb( staSt.breal_idx ) = - ones( length( staSt.breal_idx ), 1 );
staSt.ub( staSt.breal_idx ) = ones( length( staSt.breal_idx ), 1 );
staSt.scVec( staSt.breal_idx ) = RFMax / sqrt( 2 ) * ones( length( staSt.breal_idx ), 1 );

staSt.lb( staSt.bimag_idx ) = - ones( length( staSt.bimag_idx ), 1 );
staSt.ub( staSt.bimag_idx ) = ones( length( staSt.bimag_idx ), 1 );
staSt.scVec( staSt.bimag_idx ) = RFMax / sqrt( 2 ) * ones( length( staSt.bimag_idx ), 1 );

%% Specify Constraints
[ AbInitSt, AbeqInitSt, nlconIneqInitSt, nlconEqInitSt ] =...
    initializeConstraintTrackers();
binInit = true;
[ AbInitSt, AbeqInitSt, nlconIneqInitSt, nlconEqInitSt, staSt ] =...
    conAdjustSPINS(...
    opt, staSt, AbInitSt, AbeqInitSt, nlconIneqInitSt, nlconEqInitSt, binInit );

staSt = addConstraintsOpt( staSt, AbInitSt, AbeqInitSt, nlconIneqInitSt, nlconEqInitSt );

% Update possible nonlinear constraints
if ( ~isempty( opt.nlconIneqFuncs ) ) || ( ~isempty( opt.nlconEqFuncs ) )
    staSt.nonlcon = @( pSc ) nonlconCombine( pSc( : ), staSt );
else
    staSt.nonlcon = [];
end

%% Initial guess
stainit0 = zeros( staSt.numVars, 1 );
stainit0( staSt.breal_idx ) = zeros( length( staSt.breal_idx ), 1 );
stainit0( staSt.bimag_idx ) = zeros( length( staSt.bimag_idx ), 1 );
stainitSc0 = stainit0;

%% Initialize Aspatial and bspatial
sinit = struct;
sinit.kmax = initParms.kmax;
sinit.a = initParms.a;
sinit.b = initParms.b;
sinit.u = initParms.u;
sinit.v = initParms.v;
sinit.T = Tspins;
sinit.gyro = opt.gyro;

[ Aspatial, bspatial ] = generateSPINSSTAMatrices( tvec, dtvec, sinit, staSt, opt );

%% Specify Objective Function
MLSInitObjFunAnon = @( staSC ) MLSInitObjFun( staSC, Aspatial, bspatial, staSt );

%% Run initial optimization
[ stainitScopt, ~, ~, ~ ] = fmincon( MLSInitObjFunAnon, stainitSc0,...
    staSt.A, staSt.b, staSt.Aeq, staSt.beq, staSt.lb, staSt.ub, staSt.nonlcon, fminopt );
% [ stainitScopt, cost, exitflag, output ] = fmincon( MLSInitObjFunAnon, stainitSc0,...
    % staSt.A, staSt.b, staSt.Aeq, staSt.beq, staSt.lb, staSt.ub, staSt.nonlcon, fminopt );

stainitopt = staSt.scVec .* stainitScopt;

%% Prepare Second optimization
varArray = {...
        "breal", numTimePoints * opt.numXYCoils;...
        "bimag", numTimePoints * opt.numXYCoils;...
        "kmax", 1;...
        "a", 1;...
        "b", 1;...
        "u", 1;...
        "v", 1;...
        };

[ staSt, ~, ~, ~, ~, ~, ~, ~ ] =...
        processVarOrganization( staSt, varArray );

staSt.param_idx = [...
    staSt.kmax_idx;...
    staSt.a_idx;...
    staSt.b_idx;...
    staSt.u_idx;...
    staSt.v_idx;...
    ];

staSt.lb = -inf * ones( staSt.numVars, 1 );
staSt.ub = inf * ones( staSt.numVars, 1 );
scVec = zeros( staSt.numVars, 1 );

staSt.lb( staSt.breal_idx ) = - ones( length( staSt.breal_idx ), 1 );
staSt.ub( staSt.breal_idx ) = ones( length( staSt.breal_idx ), 1 );
scVec( staSt.breal_idx ) = RFMax / sqrt( 2 ) * ones( length( staSt.breal_idx ), 1 );

staSt.lb( staSt.bimag_idx ) = - ones( length( staSt.bimag_idx ), 1 );
staSt.ub( staSt.bimag_idx ) = ones( length( staSt.bimag_idx ), 1 );
scVec( staSt.bimag_idx ) = RFMax / sqrt( 2 ) * ones( length( staSt.bimag_idx ), 1 );

staSt.lb( staSt.param_idx ) = - ones( length( staSt.param_idx ), 1 );
staSt.ub( staSt.param_idx ) = ones( length( staSt.param_idx ), 1 );

scVec( staSt.kmax_idx ) = ( max_kmax - min_kmax ) / 2;
scVec( staSt.a_idx ) = ( max_a - min_a ) / 2;
scVec( staSt.b_idx ) = ( max_b - min_b ) / 2;
scVec( staSt.u_idx ) = ( max_u - min_u ) / 2;
scVec( staSt.v_idx ) = ( max_v - min_v ) / 2;

staSt.scVec = scVec;

staSt.min_kmax = min_kmax;
staSt.max_kmax = max_kmax;
staSt.min_a = min_a;
staSt.max_a = max_a;
staSt.min_b = min_b;
staSt.max_b = max_b;
staSt.min_u = min_u;
staSt.max_u = max_u;
staSt.min_v = min_v;
staSt.max_v = max_v;

%% Place results as initial guess for second optimization
sta0 = zeros( staSt.numVars, 1 );
sta0( staSt.breal_idx ) = stainitopt( staSt.breal_idx );
sta0( staSt.bimag_idx ) = stainitopt( staSt.bimag_idx );
sta0( staSt.kmax_idx ) = initParms.kmax - scVec( staSt.kmax_idx ) - min_kmax;
sta0( staSt.a_idx ) = initParms.a - scVec( staSt.a_idx ) - min_a;
sta0( staSt.b_idx ) = initParms.b - scVec( staSt.b_idx ) - min_b;
sta0( staSt.u_idx ) = initParms.u - scVec( staSt.u_idx ) - min_u;
sta0( staSt.v_idx ) = initParms.v - scVec( staSt.v_idx ) - min_v;

staSc0 = sta0 ./ scVec;

%% Specify Second Objective Function Constraints
[ AbSt, AbeqSt, nlconIneqSt, nlconEqSt ] = initializeConstraintTrackers();

binInit = false;
[ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] =...
    conAdjustSPINS( opt, staSt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt, binInit );

staSt = addConstraintsOpt( staSt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

% Update possible nonlinear constraints
if ~isempty( opt.nlconEqFuncs ) || ~isempty( opt.nlconIneqFuncs )
    staSt.nonlcon = @( pSc ) nonlconCombine( pSc( : ), staSt );
else
    staSt.nonlcon = [];
end

%% Specify Objective Function
MLSObjFunAnon = @( sta ) MLSObjFun( sta, staSt, opt );

%% Run final optimization
fminopt.MaxIterations = finalFminMaxIterations;

[ staScopt, cost, exitflag, output ] = fmincon( MLSObjFunAnon, staSc0, staSt.A, staSt.b, staSt.Aeq, staSt.beq, staSt.lb, staSt.ub, staSt.nonlcon, fminopt );
staopt = staSt.scVec .* staScopt;

output.fval = cost;
output.exitflag = exitflag;

%% Post-process sopt
sopt = struct;
sopt.kmax = min_kmax + ( staopt( staSt.kmax_idx ) + scVec( staSt.kmax_idx ) );
sopt.a = min_a + ( staopt( staSt.a_idx ) + scVec( staSt.a_idx ) );
sopt.b = min_b + ( staopt( staSt.b_idx ) + scVec( staSt.b_idx ) );
sopt.u = min_u + ( staopt( staSt.u_idx ) + scVec( staSt.u_idx ) );
sopt.v = min_v + ( staopt( staSt.v_idx ) + scVec( staSt.v_idx ) );
sopt.T = Tspins;
sopt.gyro = opt.gyro;

%% Get Gradients
grad_slew_init = ( Gxyz_fn( 0, sopt ) / initSlewTime ) * tvec( init_slew_i:init_slew_f );
grad_spins_int = Gxyz_fn( tvec( spins_int_idx ) - initSlewTime, sopt );
grad_slew_final = -( Gxyz_fn( sopt.T - finalSlewTime, sopt ) / finalSlewTime ) * ( tvec( final_slew_i:final_slew_f ) - pulseLength );

Gopt = [ grad_slew_init, grad_spins_int, grad_slew_final ];

%% Rearrange k-space vector
Kopt_noadj = kxyz_fn( tvec - initSlewTime, sopt );

% Adjust for the final slew
Kopt_adjfinalslew = Kopt_noadj;
Kopt_adjfinalslew( :, spins_int_idx ) = Kopt_adjfinalslew( :, spins_int_idx )...
    - kxyz_fn( sopt.T-finalSlewTime, sopt )...
    - ( opt.gyro * 0.5 * finalSlewTime ) * Gxyz_fn( sopt.T-finalSlewTime, sopt );
Kopt_adjfinalslew( :, final_slew_i:final_slew_f ) = ...
    - ( ( opt.gyro * Gxyz_fn( sopt.T - finalSlewTime, sopt ) ) / ( 2 * finalSlewTime ) ) *...
    ( ( tvec( final_slew_i:final_slew_f ) - pulseLength ).^2 );
Kopt_adjfinalslew = Kopt_adjfinalslew( :, spins_int_i:final_slew_f  );

% Adjust for initial slew
G0 = Gxyz_fn( 0, sopt );
deltak_initslew = ( opt.gyro / ( 2 * initSlewTime ) ) * G0 * ( tvec( init_slew_i:init_slew_f ).^2 - initSlewTime^2 );
Kopt_initslew = deltak_initslew + kxyz_fn( 0, sopt );

Kopt = [ Kopt_initslew, Kopt_adjfinalslew ];

%% Get RF waveforms
RFoptvec = complex( staopt( staSt.breal_idx ), staopt( staSt.bimag_idx ) );
RFopt = transpose( reshape( RFoptvec, [ numTimePoints, opt.numXYCoils ] ) );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ cost, gradcost ] = MLSObjFun( staSc, staSt, opt )

%% Get tvec
tvec = staSt.tvec;
dtvec = staSt.dtvec;

%% Scale parameters
sta = staSt.scVec .* staSc;

%% Get SPINS parameters
kmax = staSt.min_kmax + ( sta( staSt.kmax_idx ) + staSt.scVec( staSt.kmax_idx ) );
a = staSt.min_a + ( sta( staSt.a_idx ) + staSt.scVec( staSt.a_idx ) );
b = staSt.min_b + ( sta( staSt.b_idx ) + staSt.scVec( staSt.b_idx ) );
u = staSt.min_u + ( sta( staSt.u_idx ) + staSt.scVec( staSt.u_idx ) );
v = staSt.min_v + ( sta( staSt.v_idx ) + staSt.scVec( staSt.v_idx ) );

s.kmax = kmax;
s.a = a;
s.b = b;
s.u = u;
s.v = v;
s.T = staSt.Tspins;
s.gyro = opt.gyro;

%% Prepare Matrices
[ Aspatial, bspatial ] = generateSPINSSTAMatrices( tvec, dtvec, s, staSt, opt );

RFrealvec = sta( staSt.breal_idx );
RFimagvec = sta( staSt.bimag_idx );
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
    Xadj = ( opt.pos( :, 1 ) ) * 1j;
    Yadj = ( opt.pos( :, 2 ) ) * 1j;
    Zadj = ( opt.pos( :, 3 ) ) * 1j;

    % Reshape Aspatial for gradient with respect to k
    Aspatial_rshp = reshape( Aspatial, [ opt.numPos, staSt.numTimePoints, opt.numXYCoils ] );
    Aspatial_rshp = permute( Aspatial_rshp, [ 1, 3, 2 ] );
    RF_rshp = reshape( transpose( reshape( RFcompvec, [ staSt.numTimePoints, opt.numXYCoils ] ) ), [ opt.numXYCoils, 1, staSt.numTimePoints ] );
    Aspatial_RF_rshp = squeeze( pagemtimes( Aspatial_rshp, RF_rshp ) );

    dAxdkx = Xadj .* Aspatial_RF_rshp;
    dAxdky = Yadj .* Aspatial_RF_rshp;
    dAxdkz = Zadj .* Aspatial_RF_rshp;
    
    dfdk = [...
        real( dfdabsAx * dAxdkx );...
        real( dfdabsAx * dAxdky );...
        real( dfdabsAx * dAxdkz );...
        ];
    
    % Get initial k-space locations
    dk_noadj_dparams = dkxyzdparam_fn( tvec - staSt.initSlewTime, s );

    % Adjust for the final slew
    dk_adjfinalslew_dparams = dk_noadj_dparams;
    dk_adjfinalslew_dparams( :, staSt.spins_int_idx, : ) = dk_adjfinalslew_dparams( :, staSt.spins_int_idx, : )...
        - ( dkxyzdparam_fn( s.T-staSt.finalSlewTime, s )...
        + ( opt.gyro * 0.5 * staSt.finalSlewTime ) * dGxyzdparam_fn( s.T-staSt.finalSlewTime, s ) );
    dk_adjfinalslew_dparams( :, staSt.final_slew_i:staSt.final_slew_f, : ) = ...
        - pagemtimes( ( ( ( opt.gyro ) / ( 2 * staSt.finalSlewTime ) ) * dGxyzdparam_fn( s.T - staSt.finalSlewTime, s ) ),...
        repmat( ( tvec( staSt.final_slew_i:staSt.final_slew_f ) - s.T ).^2, [ 1, 1, length( staSt.param_idx ) ] ) );
    dk_adjfinalslew_dparams = dk_adjfinalslew_dparams( :, staSt.spins_int_i:staSt.final_slew_f, :  );

    % Adjust for initial slew
    dG0dparams = dGxyzdparam_fn( 0, s );
    dddeltak_initslew_dparams = ( opt.gyro / ( 2 * staSt.initSlewTime ) ) *...
        pagemtimes( dG0dparams, repmat( tvec( staSt.init_slew_i:staSt.init_slew_f ).^2 - staSt.initSlewTime^2, [ 1, 1, length( staSt.param_idx ) ] ) );
    dk_initslew_dparams = dddeltak_initslew_dparams + dkxyzdparam_fn( 0, s );

    dkdparams = cat( 2, dk_initslew_dparams, dk_adjfinalslew_dparams );

    dfdparams = transpose( transpose( reshape( dfdk, [ 3*staSt.numTimePoints, 1 ] ) ) *  reshape( dkdparams, [ 3*staSt.numTimePoints, length( staSt.param_idx ) ] ) );

    gradcost_unsc = zeros( size( sta, 1 ), 1 );
    gradcost_unsc( staSt.breal_idx ) = dfdRFan_real .* staSt.scVec( staSt.breal_idx );
    gradcost_unsc( staSt.bimag_idx ) = dfdRFan_imag .* staSt.scVec( staSt.bimag_idx );
    gradcost_unsc( staSt.param_idx ) = dfdparams .* staSt.scVec( staSt.param_idx );
    gradcost = gradcost_unsc / norm2bspatial;

end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ cost, gradcost ] = MLSInitObjFun( staSc, Aspatial, bspatial, staSt ) 

RFrealvec = staSt.scVec( staSt.breal_idx ) .* staSc( staSt.breal_idx );
RFimagvec = staSt.scVec( staSt.bimag_idx ) .* staSc( staSt.bimag_idx );
RFcompvec = complex( RFrealvec, RFimagvec );

Ax = Aspatial * RFcompvec;
absAx = abs( Ax );
magb = abs( bspatial );
absAxmb = absAx - magb;

norm2resid = norm( absAxmb, 2 )^2;
norm2bspatial = norm( magb, 2 )^2;
cost = norm2resid / norm2bspatial;

if nargout > 1
    
    % pre-calculate values
    if any( staSc ~= 0 )
        dfdabsAx = 2 * transpose( absAxmb .* ( conj( Ax ) ./ absAx ) );
    else
        dfdabsAx = 2 * transpose( absAxmb );
    end
    
    % deal with RF 
    dfdRFan_comp = transpose( conj( dfdabsAx * Aspatial ) );
    dfdRFan_real = real( dfdRFan_comp );
    dfdRFan_imag = imag( dfdRFan_comp );

    gradcost_unsc = zeros( size( staSc, 1 ), 1 );
    gradcost_unsc( staSt.breal_idx ) = dfdRFan_real .* staSt.scVec( staSt.breal_idx );
    gradcost_unsc( staSt.bimag_idx ) = dfdRFan_imag .* staSt.scVec( staSt.bimag_idx );

    gradcost = gradcost_unsc / norm2bspatial;
end

end
% ----------------------------------------------------------------------- %
