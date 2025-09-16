function [ oc, pulse, opt ] = processSPINSPulse_base( oc, pulse, opt )
% This function will initialize an SPINS base pulse and generate
% the needed constraints

%% Add parameters to opt struct
opt.gyro = 267.5e6;

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Get Initial Timing Parameters
[ opt, oc, tvec, ~, numTimePoints, pulseLength ] = processTimingFixedStep( oc, pulse, opt );
opt = generateRotFrameFrequencyVector( opt );

%% Check if constraints are present. If not, add in the needed ones
constraintList = keys( pulse.constraints );

if ~any( contains( constraintList, "RF-max", 'ignorecase', true ) )
    pulse.constraints( "RF-max" ) = 550; % V
end
% if ~any( contains( constraintList, "RF-slew-rate", 'ignorecase', true ) )
%     pulse.constraints( "RF-slew-rate" ) = 3e7; % V/s
% end
if ~any( contains( constraintList, "grad-max", 'ignorecase', true ) )
    pulse.constraints( "grad-max" ) = 25e-3; % T/m
end
if ~any( contains( constraintList, "grad-slew-rate", 'ignorecase', true ) )
    pulse.constraints( "grad-slew-rate" ) = 200; % T/m/sec
end
if opt.numZCoils > 0
    if ~any( contains( constraintList, "shim-max", 'ignorecase', true ) )
        pulse.constraints( "shim-max" ) = 30; % T/m
    end
    if ~any( contains( constraintList, "shim-slew-rate", 'ignorecase', true ) )
        pulse.constraints( "shim-slew-rate" ) = 3e6; % T/m/sec
    end
end

%% Determine slew limits
numSigDigitsRound = 6;

% determine the initial slew time
if isfield( pulse, "initSlewTime" )
    initSlewTime = round( round( pulse.initSlewTime/opt.dt ) * opt.dt, numSigDigitsRound, 'significant');
else
    initSlewTime = 30e-6; % 10 us minimum slew time
end

% determine the final slew time
if isfield( pulse, "finalSlewTime" )
    finalSlewTime = round( round( pulse.finalSlewTime/opt.dt ) * opt.dt, numSigDigitsRound, 'significant');
else
    finalSlewTime = 20e-6; % 10 us minimum slew time
end

%% Add slew rate information to opt struct
if opt.numZCoils > 0
    opt.shimSlewRate = pulse.constraints( 'shim-slew-rate' );
end
opt.gradSlewRate = pulse.constraints( 'grad-slew-rate' );

%% Determine points of different periods
idx_tol = opt.dt/100;
init_slew_i = find( ( tvec ) >= ( 0 - idx_tol ) , 1, 'first' );
init_slew_f = find( ( tvec ) <= ( initSlewTime + idx_tol ), 1, 'last' );
spins_i = find( ( tvec ) >= ( initSlewTime - idx_tol ), 1, 'first' );
spins_f = find( ( tvec ) <= ( pulseLength + idx_tol ), 1, 'last' );
final_slew_i = find( ( tvec ) >= ( pulseLength - finalSlewTime - idx_tol ), 1, 'first' );
final_slew_f = find( ( tvec ) <= ( pulseLength + idx_tol ), 1, 'last' );

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
spins_int_num = (final_slew_i - 1) - spins_i + 1;
Tspins = pulseLength - initSlewTime;

%% Assign to opt struct
opt.initSlewTime = initSlewTime;
opt.init_slew_i = init_slew_i;
opt.init_slew_f = init_slew_f;
opt.init_slew_sc = init_slew_sc;
opt.init_slew_num = init_slew_num;

opt.finalSlewTime = finalSlewTime;
opt.final_slew_i = final_slew_i;
opt.final_slew_f = final_slew_f;
opt.final_slew_sc = final_slew_sc;
opt.final_slew_num = final_slew_num;

opt.spins_i = spins_i;
opt.spins_f = spins_f;
opt.spins_int_i = spins_int_i;
opt.spins_int_f = spins_int_f;
opt.spins_int_idx = uint32( spins_int_i:spins_int_f );
opt.spins_int_num = spins_int_num;
opt.spins_num = spins_num;

opt.Tspins = Tspins;

opt.tvec = tvec;

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "kmax", 1;...
        "a", 1;...
        "b", 1;...
        "u", 1;...
        "v", 1;...
        "breal", numTimePoints * opt.numXYCoils;...
        "bimag", numTimePoints * opt.numXYCoils;...
        "shim", numTimePoints * opt.numZCoils;...
        };

    varArray = checkShimVariables( varArray, opt );

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );
    
    opt.estMaxActiveVarsTimeStep = 2 * opt.numXYCoils + 5 + opt.numZCoils;

elseif strcmpi( opt.structtype, "val" )
    
    varNames = opt.varNames;
    % varAmts = opt.varAmts;
    % varNumDictionary = opt.varNumDictionary;
    % varCumSum = opt.varCumSum;
    varIdxs = opt.varIdxs;
    numVars = opt.numVars;

    for vv = 1:length( varNames )
        opt.( sprintf( "%s_idx", replace(varNames( vv ), "-", "_" ) ) ) = varIdxs{ vv };
    end
end

%% Process constraints and create variable scaling
scVec = ones( numVars, 1 );
lb = -inf * ones( numVars, 1 );
ub = inf * ones( numVars, 1 );

[ AbSt, AbeqSt, nlconIneqSt, nlconEqSt ] = initializeConstraintTrackers();

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
lb( opt.kmax_idx ) = -1;
ub( opt.kmax_idx ) = 1;
scVec( opt.kmax_idx ) = ( max_kmax - min_kmax ) / 2;
opt.min_kmax = min_kmax;
opt.max_kmax = max_kmax;

% u
if isfield( pulse, "min_u" )
    min_u = pulse.min_u;
else
    min_u = ( 2*pi * 1000 ); % rad/m
end
if isfield( pulse, "max_u" )
    max_u = pulse.max_u;
else
    max_u = ( 2*pi * 10 * 1000 ); % rad/m
end
lb( opt.u_idx ) = -1;
ub( opt.u_idx ) = 1;
scVec( opt.u_idx ) = ( max_u - min_u ) / 2;
opt.min_u = min_u;
opt.max_u = max_u;

% v
if isfield( pulse, "min_v" )
    min_v = pulse.min_v;
else
    min_v = ( 2*pi * 1000 ); % rad/m
end
if isfield( pulse, "max_v" )
    max_v = pulse.max_v;
else
    max_v = ( 2*pi * 10 * 1000 ); % rad/m
end
lb( opt.v_idx ) = -1;
ub( opt.v_idx ) = 1;
scVec( opt.v_idx ) = ( max_v - min_v ) / 2;
opt.min_v = min_v;
opt.max_v = max_v;

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
lb( opt.a_idx ) = -1;
ub( opt.a_idx ) = 1;
scVec( opt.a_idx ) = ( max_a - min_a ) / 2;
opt.min_a = min_a;
opt.max_a = max_a;

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
lb( opt.b_idx ) = -1;
ub( opt.b_idx ) = 1;
scVec( opt.b_idx ) = ( max_b - min_b ) / 2;
opt.min_b = min_b;
opt.max_b = max_b;

opt.param_idx = [ opt.kmax_idx; opt.a_idx; opt.b_idx; opt.u_idx; opt.v_idx ];

%% Need to handle maximum and magnitude constraints first to get scaling vector
pulse = checkShimConstraints( pulse, opt );
constraintList = keys( pulse.constraints );
magConstraints = true( size( constraintList ) );
for cc = 1:length( constraintList )
    switch constraintList( cc )
        case { 'RF-max' } 
            RFMax = pulse.constraints( 'RF-max' );
            RFMaxidx = find( contains( varNames, [ "breal"; "bimag" ],...
                'ignorecase', true ) );
            
            opt.RFMax_constr = RFMax / sqrt( 2 );

            for ii = 1:length( RFMaxidx )
                varIdx = varIdxs{ RFMaxidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = RFMax/sqrt( 2 );
            end
        case { 'shim-max' }

            shimMax = pulse.constraints( 'shim-max' );
            shimMaxidx = find( contains( varNames,...
                "shim",...
                'ignorecase', true ) );

            opt.shimMax_constr = shimMax;

            for ii = 1:length( shimMaxidx )
                varIdx = varIdxs{ shimMaxidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = shimMax;
            end
        otherwise
            % assume constraint is specified linear inequality or convex constraint and pass to
            % next stage. If error, will catch during next stage.
            magConstraints( cc ) = false;
            continue
    end
end

% Add scVec to opt struct for convex constraints
% Also define the ub and lb
opt.scVec = scVec;
opt.ub = ub;
opt.lb = lb;

%% Handle convex constraints (including inequality and equality constraints)
% set default impedance for power calculations
if ~isfield( pulse, 'Z0' )
    pulse.Z0 = 50; % default impedance is 50 Ω
    opt.Z0 = pulse.Z0;
else
    opt.Z0 = pulse.Z0;
end
% set duty cycle for power calculations
if ~isfield( pulse, 'dutyCycle' )
    pulse.dutyCycle = 1; % default duty cycle is 100%
    opt.dutyCycle = pulse.dutyCycle;
else
    opt.dutyCycle = pulse.dutyCycle;
end

convConstraintList = constraintList( ~magConstraints );

for cc = 1:length( convConstraintList )
    switch convConstraintList( cc )
        case { 'grad-max' }

            gradMax = pulse.constraints( 'grad-max' );
            opt.gradMax_constr = gradMax;
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseGradMax, 3 );

        case { 'grad-slew-rate' }

            gradSlewRate = pulse.constraints( 'grad-slew-rate' );
            opt.gradSlewRate_constr = gradSlewRate;
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseGradSlewRate, 3 * spins_int_num + 2 * 3 );

        case { 'shim-total' }

            opt.shimTotal_constr = pulse.constraints( "shim-total" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseShimTotal, numTimePoints );

        case { 'total-RF-power' }

            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseTotalRFPower, 1 );

        case { 'max-RF-power' }
            
            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseMaxRFPower, 1 );

        case { 'peak-local-SAR' }
            
            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulsePeakLocalSAR, 1 );

        case { 'peak-global-SAR' }
            
            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulsePeakGlobalSAR, 1 );

        case { 'average-local-SAR' }
            
            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }
            
            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseAvgGlobalSAR, 1 );

        case { 'peak-E-10m' }

            opt.peakE10m_constr = pulse.constraints( "peak-E-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulsePeakE10m, 1 );

        case { 'peak-H-10m' }

            opt.peakH10m_constr = pulse.constraints( "peak-H-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulsePeakH10m, 1 );

        case { 'average-E-10m' }

            opt.avgE10m_constr = pulse.constraints( "average-E-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseAvgE10m, opt.numVOPs_E10m );

        case { 'average-H-10m' }

            opt.avgH10m_constr = pulse.constraints( "average-H-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintSPINSPulseAvgH10m, opt.numVOPs_H10m );

        case { 'RF-slew-rate' }
            
            opt.RFSlewRate_constr = pulse.constraints( 'RF-slew-rate' );
            [ ARFSlew, bRFSlew ] = constraintSPINSPulseRFSlew(...
                opt, pulse.constraints( 'RF-slew-rate' ) );

            AbSt = AbAddList( AbSt, "RF-slew-rate", ARFSlew, bRFSlew );

        case { 'RF-accel' }
            
            opt.RFAccel_constr = pulse.constraints( 'RF-accel' );
            [ ARFAccel, bRFAccel ] = constraintSPINSPulseRFAccel(...
                opt, pulse.constraints( 'RF-accel' ) );

            AbSt = AbAddList( AbSt, "RF-accel", ARFAccel, bRFAccel );

        case { 'shim-slew-rate' }
            
            opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
            [ AShimSlew, bShimSlew ] = constraintSPINSPulseShimSlew(...
                opt, pulse.constraints( 'shim-slew-rate' ) );

            AbSt = AbAddList( AbSt, "shim-slew-rate", AShimSlew, bShimSlew );

        case { 'shim-accel' }
            
            opt.shimAccel_constr = pulse.constraints( 'shim-accel' );
            [ AShimAccel, bShimAccel ] = constraintSPINSPulseShimAccel(...
                opt, pulse.constraints( 'shim-accel' ) );

            AbSt = AbAddList( AbSt, "shim-accel", AShimAccel, bShimAccel );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generatePlotWaveforms = @generateSPINSPlotWaveform_base;
opt.generateWaveforms = @generateSPINSWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientSPINS_base;

if opt.useGPU
    opt.gpuArrayAdjointFunction = @gpuArrayAdjointSPINS_base;
    opt = prepareGPUArrays( opt );
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'pulseLength' ) = pulseLength;
pulseinfo( 'initSlewTime' ) = initSlewTime;
pulseinfo( 'finalSlewTime' ) = finalSlewTime;
pulseinfo( 'Tspins' ) = Tspins;

pulseinfo( 'min_kmax' ) = min_kmax;
pulseinfo( 'max_kmax' ) = max_kmax;
pulseinfo( 'min_u' ) = min_u;
pulseinfo( 'max_u' ) = max_u;
pulseinfo( 'min_v' ) = min_v;
pulseinfo( 'max_v' ) = max_v;
pulseinfo( 'min_a' ) = min_a;
pulseinfo( 'max_a' ) = max_a;
pulseinfo( 'min_b' ) = min_b;
pulseinfo( 'max_b' ) = max_b;

pulse.pulseinfo = pulseinfo;

end