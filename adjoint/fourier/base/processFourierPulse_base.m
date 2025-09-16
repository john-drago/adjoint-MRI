function [ oc, pulse, opt ] = processFourierPulse_base( oc, pulse, opt )
% This function will initialize an "optimal control" pulse whereby all the
% RF and gradient is optimizable

%% Add parameters to opt struct
opt.gyro = 267.5e6;

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Get Timing Parameters
[ opt, oc, ~, ~, numTimePoints, pulseLength ] = processTimingFixedStep( oc, pulse, opt );
opt.tdom = [ 0, pulseLength ];
opt = generateRotFrameFrequencyVector( opt );

%% Handle order of Fourier series
if isfield( pulse, 'orderFourier_RF' )
    opt.orderFourier_RF = uint32( pulse.orderFourier_RF );
elseif isfield( pulse, 'orderFourier' )
    opt.orderFourier_RF = uint32( pulse.orderFourier );
else
    error( "Did not specify order of RF Fourier series" );
end
if opt.orderFourier_RF > ( ceil( numTimePoints / 2 ) - 1 )
    warning( "Reducing order of RF Fourier series to ensure fewer than number of points." )
    opt.orderFourier_RF = ( ceil( numTimePoints / 2 ) - 1 );
end

opt.numFourier_RF = uint32( 2 * opt.orderFourier_RF + 1 );

if isfield( pulse, 'orderFourier_grad' )
    opt.orderFourier_grad = uint32( pulse.orderFourier_grad );
elseif isfield( pulse, 'orderFourier' )
    opt.orderFourier_grad = uint32( pulse.orderFourier );
else
    error( "Did not specify order of grad Fourier series" );
end
if opt.orderFourier_grad > ( ceil( numTimePoints / 2 ) - 1 )
    warning( "Reducing order of grad Fourier series to ensure fewer than number of points." )
    opt.orderFourier_grad = ( ceil( numTimePoints / 2 ) - 1 );
end

opt.numFourier_grad = uint32( 2 * opt.orderFourier_grad + 1 );

if opt.numZCoils > 0
    if isfield( pulse, 'orderFourier_shim' )
        opt.orderFourier_shim = uint32( pulse.orderFourier_shim );
    elseif isfield( pulse, 'orderFourier' )
        opt.orderFourier_shim = uint32( pulse.orderFourier );
    else
        error( "Did not specify order of shim Fourier series" );
    end
    if opt.orderFourier_shim > ( ceil( numTimePoints / 2 ) - 1 )
        warning( "Reducing order of shim Fourier series to ensure fewer than number of points." )
        opt.orderFourier_shim = ( ceil( numTimePoints / 2 ) - 1 );
    end

    if ~isfield( opt, 'orderFourier_shim' ) || ( opt.orderFourier_shim == 0 )
        opt.numFourier_shim = 0;
        opt.numZCoils = 0;
    else
        opt.numFourier_shim = uint32( 2 * opt.orderFourier_shim + 1 );
    end

else
    if isfield( pulse, 'orderFourier_shim' )
        if pulse.orderFourier_shim > 0
            warning( "Removing shim optimization variables, because detected numZCoils = 0." );
        end
    end
    opt.orderFourier_shim = 0;
    opt.numFourier_shim = 0;
end

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
opt.gradSlewRate = pulse.constraints( 'grad-slew-rate' );
if opt.numZCoils > 0
    opt.shimSlewRate = pulse.constraints( 'shim-slew-rate' );
end

%% Determine if there is an RF-bandwidth constraint
% determine fundamental frequencies
opt.omega0 = ( 2 * pi ) / opt.pulseLength;
opt.f0 = 1 / opt.pulseLength;
opt.Omega0 = ( 2 * pi ) / opt.numTimePoints;
opt.F0 = 1 / opt.numTimePoints;

% can be handled explicitly by restricing number of harmonics
if any( contains( constraintList, "RF-bandwidth", 'ignorecase', true ) )

    opt.RFBandwidth_band = ( pulse.constraints( 'RF-bandwidth' ) );

    shift = true;
    fvec = transpose( fvecDFT( 1/opt.dt, opt.numTimePoints, shift ) );
    outOfBandwidth = ( abs( fvec ) >= ( opt.RFBandwidth_band + 1e-6 ) );
    opt.adjfvec = fvec( outOfBandwidth );
    opt.adjdf = 1 / opt.pulseLength;

    opt.RFBandwidth_constr = ( pulse.constraints( 'RF-max' ) * 2.5e-5 )...
        * length( opt.adjfvec ) * opt.adjdf;
    
    orderFourier_max_RF = uint32( floor( opt.RFBandwidth_band / ( opt.f0 - 1e-6 ) ) );
    
    if opt.orderFourier_RF > orderFourier_max_RF
        
        warning( "Reducing order of RF Fourier series to: %i. Will ensure RF bandwidth constraint is satisfied.", orderFourier_max_RF );

        opt.orderFourier_RF = orderFourier_max_RF;
        opt.numFourier_RF = uint32( 2 * opt.orderFourier_RF + 1 );
    end
    
end

%% Initialize Fourier samples and get frequencies
if opt.numZCoils > 0
    opt.orderFourier_max = max( [ opt.orderFourier_RF, opt.orderFourier_grad, opt.orderFourier_shim ] );
else
    opt.orderFourier_max = max( [ opt.orderFourier_RF, opt.orderFourier_grad ] );
end

[ opt.FBRF, opt.FBRFfvec ] = getFourierRFBasisMatrix( opt.orderFourier_RF, opt.tvec, opt.omega0 );

[ opt.FBgrad, opt.FBgradfvec ] = getFourierBasisMatrix( opt.orderFourier_grad, opt.tvec, opt.omega0 );

if opt.numZCoils
    
    [ opt.FBshim, opt.FBshimfvec ] = getFourierBasisMatrix( opt.orderFourier_shim, opt.tvec, opt.omega0 );

end

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal", ( opt.numFourier_RF ) * opt.numXYCoils;...
        "bimag", ( opt.numFourier_RF ) * opt.numXYCoils;...
        "grad", opt.numFourier_grad * 3;...
        "shim", opt.numFourier_shim * opt.numZCoils;...
        };

    varArray = checkShimVariables( varArray, opt );

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );
    
    opt.estMaxActiveVarsTimeStep = opt.numXYCoils * ( 2 * opt.numFourier_RF ) ...
        + opt.numFourier_grad * 3 + opt.numFourier_shim * opt.numZCoils;

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
            RFMaxidx = find( contains( varNames, "b",...
                'ignorecase', true ) );

            opt.RFMax_constr = RFMax / sqrt(2);

            for ii = 1:length( RFMaxidx )
                varIdx = varIdxs{ RFMaxidx( ii ) };
                % lb( varIdx ) = -1;
                % ub( varIdx ) = 1;
                scVec( varIdx ) = RFMax / sqrt( 2 );
            end

            % Set value at t=0 and t=pulseLength to be zero
            breal_idx_rshp = reshape( opt.breal_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
            bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
            
            % breal
            [ Aeqbrealif, beqbrealif ] = constraintFourierWaveformInitFinal(...
                breal_idx_rshp, opt.numFourier_RF, opt );
            AbeqSt = AbAddList( AbeqSt, "breal-if", Aeqbrealif, beqbrealif );
            
            % bimag
            [ Aeqimagif, beqbimagif ] = constraintFourierWaveformInitFinal(...
                bimag_idx_rshp, opt.numFourier_RF, opt );
            AbeqSt = AbAddList( AbeqSt, "bimag-if", Aeqimagif, beqbimagif );

            % still have to implement RF-max constraints in chebyshev basis
            magConstraints( cc ) = false;

        case { 'grad-max' }

            gradMax = pulse.constraints( 'grad-max' );
            gradMaxidx = find( contains( varNames,...
                "grad",...
                'ignorecase', true ) );
            
            opt.gradMax_constr = gradMax;

            for ii = 1:length( gradMaxidx )
                varIdx = varIdxs{ gradMaxidx( ii ) };
                % lb( varIdx ) = -1;
                % ub( varIdx ) = 1;
                scVec( varIdx ) = gradMax;
            end
            
            % Set value at t=0 and t=pulseLength to be zero
            grad_idx_rshp = reshape( opt.grad_idx, [ opt.numFourier_grad, 3 ] );

            % grad
            [ Aeqgradif, beqgradif ] = constraintFourierWaveformInitFinal(...
                grad_idx_rshp( 1:(opt.orderFourier_grad+1), : ), (opt.orderFourier_grad+1), opt );
            AbeqSt = AbAddList( AbeqSt, "grad-if", Aeqgradif, beqgradif );

            magConstraints( cc ) = false;

        case { 'shim-max' }

            shimMax = pulse.constraints( 'shim-max' );
            shimMaxidx = find( contains( varNames,...
                "shim",...
                'ignorecase', true ) );

            opt.shimMax_constr = shimMax;

            for ii = 1:length( shimMaxidx )
                varIdx = varIdxs{ shimMaxidx( ii ) };
                % lb( varIdx ) = -1;
                % ub( varIdx ) = 1;
                scVec( varIdx ) = shimMax;
            end

            % Set value at t=0 and t=pulseLength to be zero
            shim_idx_rshp = reshape( opt.shim_idx, [ opt.numFourier_shim, opt.numZCoils ] );

            % shim
            [ Aeqshimif, beqshimif ] = constraintFourierWaveformInitFinal(...
                shim_idx_rshp( 1:(opt.orderFourier_shim+1), : ), (opt.orderFourier_shim+1), opt );
            AbeqSt = AbAddList( AbeqSt, "shim-if", Aeqshimif, beqshimif );

            magConstraints( cc ) = false;
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
        case { 'RF-max' }

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierRFMax, 2*opt.numXYCoils );

        case { 'RF-slew-rate' }
            
            opt.RFSlewRate_constr = pulse.constraints( 'RF-slew-rate' ); 
            opt.D_RF = getFourierRFBasisDerivativeMatrix( opt.orderFourier_RF, opt.FBRF, opt.omega0 );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierRFSlewRate, 2*opt.numXYCoils );

        case { 'RF-accel' }

            opt.RFAccel_constr = pulse.constraints( 'RF-accel' ); 
            opt.D_RF = getFourierRFBasisDerivativeMatrix( opt.orderFourier_RF, opt.FBRF, opt.omega0 );
            opt.D2_RF = getFourierRFBasisDerivativeMatrix( opt.orderFourier_RF, opt.D_RF, opt.omega0 );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierRFAccel, 2*opt.numXYCoils );
        
        case { 'grad-max' }
            
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierGradMax, 3 );

        case { 'grad-slew-rate' }

            opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );
            opt.D_grad = getFourierBasisDerivativeMatrix( opt.orderFourier_grad, opt.tvec, opt.omega0 );
            maxExtrema = 10;
            maxExtrema = min( [ maxExtrema, opt.orderFourier_grad - 1 ] );
            maxExtrema = max( [ 1, maxExtrema ] );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @( pSc, opt ) constraintFourierGradSlewRate( pSc, opt, maxExtrema ), 3 * maxExtrema );

        case { 'grad-accel' }

            opt.gradAccel_constr = pulse.constraints( 'grad-accel' );
            opt.D2_grad = getFourierBasis2DerivativeMatrix( opt.orderFourier_grad, opt.tvec, opt.omega0 );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierGradAccel, 3 );

        case { 'shim-max' }

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierShimMax, opt.numZCoils );

        case { 'shim-slew-rate' }
        
            opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' ); 
            opt.D_shim = getFourierBasisDerivativeMatrix( opt.orderFourier_shim, opt.tvec, opt.omega0 );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierShimSlewRate, opt.numZCoils );

        case { 'shim-accel' }
        
            opt.shimAccel_constr = pulse.constraints( 'shim-accel' );
            opt.D2_shim = getFourierBasis2DerivativeMatrix( opt.orderFourier_shim, opt.tvec, opt.omega0 );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierShimAccel, opt.numZCoils );

        case { 'shim-total' }
        
            opt.shimTotal_constr = pulse.constraints( 'shim-total' ); 
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierShimTotal, 1 );

        case { 'total-RF-power' }

            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierMaxRFPower, 1 );

        case { 'RF-bandwidth' }

        case { 'peak-local-SAR' }

            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierPeakLocalSAR, 1 );

        case { 'peak-global-SAR' }

            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierPeakGlobalSAR, 1 );

        case { 'average-local-SAR' }

            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintFourierAvgGlobalSAR, 1 );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generateWaveforms = @generateFourierWaveform_base;
opt.generatePlotWaveforms = @generateFourierPlotWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientFourier_base;

if matches( opt.structtype, 'opt' )
    if opt.useGPU
        opt.gpuArrayAdjointFunction = @gpuArrayAdjointFourier_base;
        opt = prepareGPUArrays( opt );
        opt.gpu_FBRF = gpuArray( opt.FBRF );
        opt.gpu_FBgrad = gpuArray( opt.FBgrad );
        if opt.numZCoils > 0
            opt.gpu_FBshim = gpuArray( opt.FBshim );
        end
    end
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'numTimePoints' ) = numTimePoints;
pulseinfo( 'dt' ) = opt.dt;
pulseinfo( 'dutyCycle' ) = opt.dutyCycle;
pulseinfo( 'orderFourier_max' ) = opt.orderFourier_max;
pulseinfo( 'orderFourier_RF' ) = opt.orderFourier_RF;
pulseinfo( 'orderFourier_grad' ) = opt.orderFourier_grad;
if opt.numZCoils > 0
    pulseinfo( 'orderFourier_shim' ) = opt.orderFourier_shim;
else
    pulseinfo( 'orderFourier_shim' ) = 0;
end

pulse.pulseinfo = pulseinfo;

end