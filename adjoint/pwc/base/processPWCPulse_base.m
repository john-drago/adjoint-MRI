function [ oc, pulse, opt ] = processPWCPulse_base( oc, pulse, opt )
% This function will initialize an "optimal control" pulse whereby all the
% RF and gradient is optimizable

%% Add parameters to opt struct
opt.gyro = 267.5e6;

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Get Timing Parameters
[ opt, oc, ~, ~, numTimePoints ] = processTimingFixedStep( oc, pulse, opt );
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
        pulse.constraints( "shim-max" ) = 30; % Amp/turn
    end
    if ~any( contains( constraintList, "shim-slew-rate", 'ignorecase', true ) )
        pulse.constraints( "shim-slew-rate" ) = 3e6; % Amp/turn/sec
    end
end

%% Add slew rate information to opt struct
opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );
if opt.numZCoils > 0
    opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
end

%% Assign to opt struct

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal", numTimePoints * opt.numXYCoils;...
        "bimag", numTimePoints * opt.numXYCoils;...
        "grad",  numTimePoints * 3;...
        "shim",  numTimePoints * opt.numZCoils;...
        };

    varArray = checkShimVariables( varArray, opt );

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );

    opt.estMaxActiveVarsTimeStep = 2 * opt.numXYCoils + 3 + opt.numZCoils;

elseif strcmpi( opt.structtype, "val" )
    
    varNames = opt.varNames;
    % varAmts = opt.varAmts;
    % varNumDictionary = opt.varNumDictionary;
    % varCumSum = opt.varCumSum;
    varIdxs = opt.varIdxs;
    numVars = opt.numVars;
end

%% Assign indices to opt struct for faster indexing in functions
for vv = 1:length( varNames )
    opt.( sprintf( "%s_idx", replace(varNames( vv ), "-", "_" ) ) ) = varIdxs{ vv };
end

%% Process constraints and create variable scaling
scVec = ones( numVars, 1 );
lb = -inf * ones( numVars, 1 );
ub = inf * ones( numVars, 1 );

[ AbSt, AbeqSt, nlconIneqSt, nlconEqSt ] = initializeConstraintTrackers();

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

        case { 'grad-max' }
            gradMax = pulse.constraints( 'grad-max' );
            gradMaxidx = find( contains( varNames,...
                "grad",...
                'ignorecase', true ) );
            
            opt.gradMax_constr = gradMax;

            for ii = 1:length( gradMaxidx )
                varIdx = varIdxs{ gradMaxidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = gradMax;
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
        case { 'total-RF-power' }

            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseMaxRFPower, 1 );

        case { 'peak-local-SAR' }

            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulsePeakLocalSAR, 1 );

        case { 'peak-global-SAR' }

            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulsePeakGlobalSAR, 1 );

        case { 'average-local-SAR' }

            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseAvgGlobalSAR, 1 );

        case { 'peak-E-10m' }

            opt.peakE10m_constr = pulse.constraints( "peak-E-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulsePeakE10m, 1 );

        case { 'peak-H-10m' }

            opt.peakH10m_constr = pulse.constraints( "peak-H-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulsePeakH10m, 1 );

        case { 'average-E-10m' }

            opt.avgE10m_constr = pulse.constraints( "average-E-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseAvgE10m, opt.numVOPs_E10m );

        case { 'average-H-10m' }

            opt.avgH10m_constr = pulse.constraints( "average-H-10m" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseAvgH10m, opt.numVOPs_H10m );
        
        case { 'RF-bandwidth' }

            opt.RFBandwidth_band = ( pulse.constraints( 'RF-bandwidth' ) );

            % determine DFT matrix
            DFTmtx = fftshift( fft( eye( opt.numTimePoints ), [], 1 ), 1 );
            shift = true;
            fvec = transpose( fvecDFT( 1/opt.dt, opt.numTimePoints, shift ) );
            outOfBandwidth = ( abs( fvec ) >= ( opt.RFBandwidth_band + 1e-6 ) );
            opt.adjDFTmtx = DFTmtx( outOfBandwidth, : ) / opt.numTimePoints;
            opt.adjfvec = fvec( outOfBandwidth );
            opt.adjdf = 1 / opt.pulseLength;

            opt.RFBandwidth_constr = ( RFMax * 1.0e-5 ) * length( opt.adjfvec ) * opt.adjdf;

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseRFBandwidth, opt.numXYCoils );

        case { 'RF-slew-rate' }

            opt.RFSlewRate_constr = pulse.constraints( 'RF-slew-rate' );
            [ ARFSlew, bRFSlew ] = constraintPWCPulseRFSlew(...
                opt, pulse.constraints( 'RF-slew-rate' ) );

            AbSt = AbAddList( AbSt, "RF-slew-rate", ARFSlew, bRFSlew );

        case { 'RF-accel' }

            opt.RFAccel_constr = pulse.constraints( 'RF-accel' );
            [ ARFAccel, bRFAccel ] = constraintPWCPulseRFAccel(...
                opt, pulse.constraints( 'RF-accel' ) );

            AbSt = AbAddList( AbSt, "RF-accel", ARFAccel, bRFAccel );

        case { 'grad-slew-rate' }

            opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );
            [ AGradSlew, bGradSlew ] = constraintPWCPulseGradSlew(...
                opt, pulse.constraints( 'grad-slew-rate' ) );

            AbSt = AbAddList( AbSt, "grad-slew-rate", AGradSlew, bGradSlew );

        case { 'grad-accel' }

            opt.gradAccel_constr = pulse.constraints( 'grad-accel' );
            [ AGradAccel, bGradAccel ] = constraintPWCPulseGradAccel(...
                opt, pulse.constraints( 'grad-accel' ) );

            AbSt = AbAddList( AbSt, "grad-accel", AGradAccel, bGradAccel );

        case { 'shim-slew-rate' }

            opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
            [ AShimSlew, bShimSlew ] = constraintPWCPulseShimSlew(...
                opt, pulse.constraints( 'shim-slew-rate' ) );

            AbSt = AbAddList( AbSt, "shim-slew-rate", AShimSlew, bShimSlew );

        case { 'shim-accel' }

            opt.shimAccel_constr = pulse.constraints( 'shim-accel' );
            [ AShimAccel, bShimAccel ] = constraintPWCPulseShimAccel(...
                opt, pulse.constraints( 'shim-accel' ) );

            AbSt = AbAddList( AbSt, "shim-accel", AShimAccel, bShimAccel );

        case { 'shim-total' }

            opt.shimTotal_constr = pulse.constraints( "shim-total" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPWCPulseShimTotal, numTimePoints );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generateWaveforms = @generatePWCWaveform_base;
opt.generatePlotWaveforms = @generatePWCPlotWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientPWC_base;

if opt.useGPU
    opt.gpuArrayAdjointFunction = @gpuArrayAdjointPWC_base;
    opt = prepareGPUArrays( opt );
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'numTimePoints' ) = opt.numTimePoints;
pulseinfo( 'dt' ) = opt.dt;

pulse.pulseinfo = pulseinfo;

end