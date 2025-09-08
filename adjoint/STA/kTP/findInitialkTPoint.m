function [ RF, G, K, AE, bE ] = findInitialkTPoint(...
    opt, pulse, RFMax, gradSlewRate, p, maxIterations )

if nargin < 3
    RFMax = 550;
end
if nargin < 4
    gradSlewRate = 200;
end
if nargin < 5
    p = 8; % set p hyperparameter that determines how many potential kT-points to evaluate
end
if nargin < 6
    maxIterations = 500;
end

%% Determine if the kTP parameters are valid
num_kTP = pulse.num_kTP;
dt_tol = 1e-6;

if num_kTP < 1
    pulse.num_kTP = 1;
    num_kTP = pulse.num_kTP;
    warning( "Changed number of kT-points to be 1." )
end

estkTPlength = num_kTP * pulse.RFLength + ( num_kTP - 1 ) * pulse.blipLength;
if abs( estkTPlength - pulse.length ) / pulse.length > dt_tol
    warning( "Number of kT-points with RF and blip length are not compatible for specified pulse length (pulse.length)." );
    warning( "Changing pulse length to %g ms", estkTPlength*1e3 );
    pulse.length = estkTPlength;
end

%% Scale RF constraints
opt.scVec( opt.breal_idx ) = RFMax / sqrt( 2 );
opt.scVec( opt.bimag_idx ) = RFMax / sqrt( 2 );

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
FOVz = 25e-2; % 20 cm FOV in z

dkx = 1 / FOVx; % units 1/m
dky = 1 / FOVy; % units 1/m
dkz = 1 / FOVz; % units 1/m

Nx = 64;
Ny = 64;
Nz = 64;

% Potential kT-points to sample
kx = ( ceil( -Nx/2 ) + ( 1:Nx ) ) * dkx;
ky = ( ceil( -Ny/2 ) + ( 1:Ny ) ) * dky;
kz = ( ceil( -Nz/2 ) + ( 1:Nz ) ) * dkz;

ikx0 = find( kx == 0 );
iky0 = find( ky == 0 );
ikz0 = find( kz == 0 );

[ Kx, Ky, Kz ] = ndgrid( kx, ky, kz );
% Kr = sqrt( Kx.^2 + Ky.^2 + Kz.^2 );

%% Determine indices for the different kT-points to sample
domSize = size( Kx );
[ I, J, K ] = ndgrid(...
    1:domSize( 1 ),...
    1:domSize( 2 ),...
    1:domSize( 3 )...
    );

ind_3D = zeros( domSize );
ind_3D( : ) = 1:prod( domSize );

%% Calculate universal offset points
dk_tol = max( [ dkx, dky, dkz ] ) / 8;
dkmax = abs( opt.gyro/(2*pi) * ( pulse.blipLength / 2 )^2 * gradSlewRate ) - dk_tol; % units 1/m
dkxmax = floor(dkmax/dkx) * dkx;
dkymax = floor(dkmax/dky) * dky;
dkzmax = floor(dkmax/dkz) * dkz;

kxo = -dkxmax : dkx : dkxmax;
kyo = -dkymax : dky : dkymax;
kzo = -dkzmax : dkz : dkzmax;

[ Kxo, Kyo, Kzo ] = ndgrid( kxo, kyo, kzo );
rko = [ Kxo( : ), Kyo( : ), Kzo( : ) ];
norm_rko = vecnorm( rko, 'inf', 2 );

[ Io, Jo, Ko ] = ndgrid( 1:length(kxo), 1:length(kyo), 1:length(kzo) );

% Sort based on magnitude of distance from point of interest
[ ~, sort_norm_rko ] = sort( norm_rko, 1, 'ascend' );

epsround = eps( 1e1 );
ijkoi = [ find( abs( kxo ) < epsround ), find( abs( kyo ) < epsround ), find( abs( kzo ) < epsround) ]; % find points that are at ro = [ 0, 0, 0 ]
ijko = [ Io( : ), Jo( : ), Ko( : ) ];
ijko =  ijkoi - ijko( sort_norm_rko, : );

ijk = [ I( : ), J( : ), K( : ) ];

%% Run optimal matching pursuit (OMP) algorithm
RFLength = pulse.RFLength;
% minRFSlewTime = pulse.minRFSlewTime;
blipLength = pulse.blipLength;
pulseLength = pulse.length;

% initialize optimization parameters for solve
fminopt = optimoptions( "fmincon" );
fminopt.ConstraintTolerance = 1e-9;
fminopt.FunctionTolerance = 1e-8;
fminopt.Display = "off";
fminopt.SpecifyConstraintGradient = true;
fminopt.SpecifyObjectiveGradient = true;
fminopt.StepTolerance = 1e-14;
fminopt.MaxFunctionEvaluations = inf;

fminopt.Algorithm = "active-set";
fminopt.RelLineSrchBnd = 1e-2;
fminopt.RelLineSrchBndDuration = inf;
fminopt.TolConSQP = 1e-10;
fminopt.MaxIterations = maxIterations;

% initialization
kk = 1;
K = [ 0; 0; 0 ]; % start with K = [0; 0; 0]; as the initial kT-point
ijklast = [ ikx0, iky0, ikz0 ];

% determine initial RF
curr_num_kTP = 1;
[ AE, bE ] = generatekTPSTAMatrices( K, pulse, opt );

RF0p = [];
[ RF, resid ] = solveOMPMatriceskTP( AE, bE, RF0p, opt, curr_num_kTP, fminopt );
RF0p = RF;

% update target magnetization
bE = abs( bE ) .* exp( 1j * angle( AE*RF ) );

% Run Loop
while kk < num_kTP

    %% Trying to find next kTP
    kk = kk + 1;
    
    %% Determine Points to evaluate based on gradient limits
    % Find latest kT-point
    % klast = K( :, 1 );
    kTTime = pulseLength...
        - ( RFLength/2 )...
        - ( kk - 1 ) * ( RFLength ) - ( kk - 1 ) * blipLength;

    % Determine ball of shortest norm around point of interest
    ijk_sort_ball = ijklast + ijko;

    % Ensure that indices are within domain bounds
    ijk_sort_lb = all( ijk_sort_ball >= 1, 2 );
    ijk_sort_ub = all(...
        [...
        ( ijk_sort_ball(:, 1) <= domSize( 1 ) ),...
        ( ijk_sort_ball(:, 2) <= domSize( 2 ) ),...
        ( ijk_sort_ball(:, 3) <= domSize( 3 ) ) ], 2 );
    ijk_sort_valid = ijk_sort_lb & ijk_sort_ub;
    ijk_sort_ball = ijk_sort_ball( ijk_sort_valid, : );

    lini_sort_ball = ...
        ijk_sort_ball( :, 1 ) +...
        ( ijk_sort_ball( :, 2 ) - 1 ) * domSize( 1 ) +...
        ( ijk_sort_ball( :, 3 ) - 1 ) * domSize( 1 ) * domSize( 2 );

    %% Evaluate All Potential Points for ability to reduce residual
    cullIdxLin = cullSumSquaresCorrelationkTP(...
        p, resid, b1pnormalized, kTTime, Kx, Ky, Kz, lini_sort_ball, opt, pulse );
    cullIdx = ind_3D( cullIdxLin );
    kT_trial = transpose( [ Kx( cullIdx ), Ky( cullIdx ), Kz( cullIdx ) ] );

    %% Select best point to reduce residual
    [ K, AE, bE, RF, resid, bestIdx ] =...
        evaluatePotentialkTPointsOMP( AE, bE, RF0p, K, kTTime, kT_trial, opt, pulse, fminopt );
    RF0p = RF;
    
    %% Update for next iteration
    ijklast = ijk( cullIdx( bestIdx ), : );
    
end

%% Calculate Gradient Blips
G = calculateGradientBlipskTPSTA( K, opt, pulse );

%% Calculate RF
RF = reshape( RF, [ opt.numXYCoils, num_kTP ] );

end