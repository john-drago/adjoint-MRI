function [ RF, K, G, spokes, AE, bE ] = findInitialSpokes( opt, pulse,...
    gradMax, gradSlewRate, p, roundTime, maxIterations )

if nargin < 3
    if isfield( opt, 'gradMax_constr' )
        gradMax = min( opt.gradMax_constr, 0.90 * 50e-3 );
    else
        gradMax = 0.90 * 50e-3;
    end
end

if nargin < 4
    if isfield( opt, 'gradSlewRate_constr' )
        gradSlewRate = min( opt.gradSlewRate_constr, 0.90 * 200 );
    else
        gradSlewRate = 0.90 * 200;
    end
end

if nargin < 5
    p = 8; % set p hyperparameter that determines how many potential spokes locations to evaluate
end

if nargin < 6
    roundTime = 5e-6;
end

if nargin < 7
    maxIterations = 3e2;
end

%% Get slice thickness
if ~isfield( pulse, 'sliceThickness' )
    error( "Unknown slice thickness." );
end

%% Determine num spokes
if isfield( pulse, 'numSpokes' )
    numSpokes = pulse.numSpokes;
else
    error( "Undefined number of spokes." )
end

%% Determine slice location
if ~isfield( pulse, 'sliceLocation' )
    error( "Unknown slice location." )
end

if isfield( pulse, 'sliceDirection' )
    pulse.sliceDirection = pulse.sliceDirection / norm( pulse.sliceDirection, 2 );
    if numel( pulse.sliceDirection ) ~= 3
        error( "Unsure how to process slice direction" );
    end
else
    warning( "Undefined slice direction. Setting slize direction to [0; 0; 1]" );
    pulse.sliceDirection = [ 0; 0; 1 ];
end

%% Determine spokes trajectory 
% round to nearest time
pulse.gyro = opt.gyro;
[ spokesTiming ] = determineSpokesDuration( pulse, gradMax, gradSlewRate, roundTime );

%% Determine spokes timing
spokesTiming.dt = 2.5e-6; % set integration time for spokes pulses. Midpoint rule
spokes = determineSpokesTimings( spokesTiming );

%% Scale RF constraints and get correct indices
opt.scVec = ( 1 * ( opt.RFMax_constr / sqrt(2) ) ) * ones( 2 * numSpokes * opt.numXYCoils, 1 );
opt.breal_idx = uint32( transpose( 1:( numSpokes * opt.numXYCoils ) ) );
opt.bimag_idx = uint32( ( numSpokes * opt.numXYCoils ) + transpose( 1:( numSpokes * opt.numXYCoils ) ) );

%% Normalize sensitivity matrices
b1p = complex( opt.b1preal, opt.b1pimag );
b1pmag = sqrt( sum( opt.b1preal.^2 + opt.b1pimag.^2, 1 ) );
b1pnormalized = b1p ./ b1pmag;

%% Adjust initial phase target to be that of the birdcage coil
CPangles = linspace( 0, 2*pi, (8+1) );
CPanglesComp = exp( 1j * CPangles( 1:(end-1) ) );
CPanglesFOV = repmat( CPanglesComp, [ opt.numPos, 1 ] );
CPmode = sum( b1p .* CPanglesFOV, 2 );
CPphase = angle( CPmode );
Mxytarg = sqrt( opt.Mtarg( :, 1 ).^2 + opt.Mtarg( :, 2 ).^2 ) .* exp( 1j * CPphase );
opt.Mtarg( :, 1 ) = real( Mxytarg );
opt.Mtarg( :, 2 ) = imag( Mxytarg );

%% Determine points to evaluate
FOVx = 25e-2; % 20 cm FOV in x
FOVy = 25e-2; % 20 cm FOV in y

dkx = 1 / FOVx; % units 1/m
dky = 1 / FOVy; % units 1/m

Nx = 64;
Ny = 64;

% Potential kT-points to sample
kx = ( ceil( -Nx/2 ) + ( 1:Nx ) ) * dkx;
ky = ( ceil( -Ny/2 ) + ( 1:Ny ) ) * dky;

ikx0 = find( kx == 0 );
iky0 = find( ky == 0 );

[ Kx, Ky ] = ndgrid( kx, ky );

domSize = size( Kx );

ind_2D = zeros( domSize );
ind_2D( : ) = 1:prod( domSize );

%% Calculate universal offset points for central k-spoke
if numSpokes > 1
    [ ij_beforeCenter, ijo_beforeCenter ] = getOffsetIdxs( dkx, dky, spokes.lastBlipLength, gradSlewRate, domSize, opt );
    [ ij_beforeNonCenter, ijo_beforeNonCenter ] = getOffsetIdxs( dkx, dky, spokes.spokeBlipLength, gradSlewRate, domSize, opt );
end

%% Get fmincon options for optimization
% initialize optimization parameters for solve
fminopt = optimoptions( "fmincon" );
fminopt.SpecifyConstraintGradient = true;
fminopt.SpecifyObjectiveGradient = true;
fminopt.StepTolerance = 1e-5;
fminopt.FunctionTolerance = 1e-6;
fminopt.ConstraintTolerance = 1e-14;
fminopt.OptimalityTolerance = 1e-7;
fminopt.MaxFunctionEvaluations = inf;

fminopt.MaxIterations = maxIterations;
fminopt.Display = "off";

% fminopt.Algorithm = "sqp";

fminopt.Algorithm = "active-set";
fminopt.RelLineSrchBnd = 1e-2;
fminopt.RelLineSrchBndDuration = inf;
fminopt.TolConSQP = 1e-10;

%% Run orthogonal matching pursuit (OMP) algorithm
% initialization
kk = 1;
K = [ 0; 0 ]; % start with K = [ 0; 0 ]; as the initial spoke location
ijlast = [ ikx0, iky0 ];

% determine initial RF
curr_numSpokes = 1;
[ AE, bE ] = generateSpokesSTAMatrices( K, spokes, opt );

RF0p = [];
[ RF, resid ] = solveOMPMatricesSpokes( AE, bE, RF0p, spokes, opt, curr_numSpokes, fminopt );
RF0p = RF;

% update target magnetization
bE = abs( bE ) .* exp( 1j * angle( AE*RF ) );

% Run Loop
while kk < numSpokes

    %% Trying to find next spoke idx
    kk = kk + 1;
    spokeIdx = numSpokes - kk + 1;
    
    %% Determine Points to evaluate based on gradient limits
    if spokeIdx == numSpokes
        ijko = ijo_beforeCenter;
        ij = ij_beforeCenter;
    else
        ijko = ijo_beforeNonCenter;
        ij = ij_beforeNonCenter;
    end

    % Determine ball of shortest norm around point of interest
    ij_sort_ball = ijlast + ijko;

    % Ensure that indices are within domain bounds
    ij_sort_lb = all( ij_sort_ball >= 1, 2 );
    ij_sort_ub = all(...
        [...
        ( ij_sort_ball(:, 1) <= domSize( 1 ) ),...
        ( ij_sort_ball(:, 2) <= domSize( 2 ) ) ], 2 );
    ij_sort_valid = ij_sort_lb & ij_sort_ub;
    ij_sort_ball = ij_sort_ball( ij_sort_valid, : );

    lini_sort_ball = ...
        ij_sort_ball( :, 1 ) +...
        ( ij_sort_ball( :, 2 ) - 1 ) * domSize( 1 );

    %% Evaluate All Potential Points for ability to reduce residual
    cullIdxLin = cullSumSquaresCorrelationSpokes(...
        p, resid, b1pnormalized, spokeIdx, Kx, Ky, lini_sort_ball, opt, spokes );
    cullIdx = ind_2D( cullIdxLin );
    spoke_trial = transpose( [ Kx( cullIdx ), Ky( cullIdx ) ] );

    %% Select best point to reduce residual
    [ K, AE, bE, RF, resid, bestIdx ] =...
        evaluatePotentialSpokesOMP( AE, bE, RF0p, K, spokeIdx, spoke_trial, opt, spokes,...
        fminopt);
    RF0p = RF;

    %% Update for next iteration
    ijlast = ij( cullIdx( bestIdx ), : );
    
end

%% Calculate Gradient Blips
G = calculateGradientBlipsSpokesSTA( K, opt, spokes );

%% Calculate RF
RF = reshape( RF, [ opt.numXYCoils, numSpokes ] );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ ij, ijo ] = getOffsetIdxs( dkx, dky, blipLength, gradSlewRate, domSize, opt )

[ I, J ] = ndgrid(...
    1:domSize( 1 ),...
    1:domSize( 2 )...
    );

dk_tol = max( [ dkx, dky ] ) / 8;
dkmax = abs( opt.gyro/(2*pi) * ( blipLength / 2 )^2 * gradSlewRate ) - dk_tol; % units 1/m
dkxmax = floor(dkmax/dkx) * dkx;
dkymax = floor(dkmax/dky) * dky;

kxo = -dkxmax : dkx : dkxmax;
kyo = -dkymax : dky : dkymax;

[ Kxo, Kyo ] = ndgrid( kxo, kyo );
rko = [ Kxo( : ), Kyo( : ) ];
norm_rko = vecnorm( rko, 'inf', 2 );

[ Io, Jo ] = ndgrid( 1:length(kxo), 1:length(kyo) );

% Sort based on magnitude of distance from point of interest
[ ~, sort_norm_rko ] = sort( norm_rko, 1, 'ascend' );

epsround = eps( 1e1 );
ijoi = [ find( abs( kxo ) < epsround ), find( abs( kyo ) < epsround ) ]; % find points that are at ro = [ 0, 0 ]
ijo = [ Io( : ), Jo( : ) ];
ijo =  ijoi - ijo( sort_norm_rko, : );

ij = [ I( : ), J( : ) ];

end
% ----------------------------------------------------------------------- %