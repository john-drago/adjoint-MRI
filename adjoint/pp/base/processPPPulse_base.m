function [ oc, pulse, opt ] = processPPPulse_base( oc, pulse, opt )
% This function will initialize an "optimal control" pulse whereby all the
% RF and gradient is optimizable

%% Add parameters to opt struct
opt.gyro = 267.5e6;

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Determine Timing
% First need to determine if the different timing structs are present
% RF
if isfield( pulse, 'timing_RF' )
    timing_RF = pulse.timing_RF;
elseif isfield( pulse, 'timing_PP' )
    timing_RF = pulse.timing_PP;
else
    error( "Unknown timing for RF piecewise polynomial." )
end
PPStartStop_RF = processStartStopPPTiming( timing_RF );
RF_end = PPStartStop_RF( end, 2 );

% Grad
if isfield( pulse, 'timing_grad' )
    timing_grad = pulse.timing_grad;
elseif isfield( pulse, 'timing_PP' )
    timing_grad = pulse.timing_PP;
else
    error( "Unknown timing for Grad piecewise polynomial." )
end
PPStartStop_grad = processStartStopPPTiming( timing_grad );
Grad_end = PPStartStop_grad( end, 2 );

% Shim
if opt.numZCoils > 0
    if isfield( pulse, 'timing_shim' )
        timing_shim = pulse.timing_shim;
    elseif isfield( pulse, 'timing_PP' )
        timing_shim = pulse.timing_PP;
    else
        error( "Unknown timing for Shim piecewise polynomial." )
    end
    PPStartStop_shim = processStartStopPPTiming( timing_shim );
    Shim_end = PPStartStop_shim( end, 2 );
end

% Check to make sure the timing structs end with pulse.length
timeResDigits = 6; % round to microseconds
pulse.length = round( pulse.length, timeResDigits); 
dt_tol = 1e-6;
idxtol = opt.dt/10;

if ( abs( RF_end - pulse.length ) < dt_tol ) || ( abs( Grad_end - pulse.length ) < dt_tol ) || ...
    ( opt.numZCoils > 0 && ( abs( Shim_end - pulse.length ) < dt_tol ) )
    opt.pulseLength = pulse.length;
else
    error( "Improper piecewise polynomial timing. The end of the timing vectors must coincide with pulse.length" );
end

%% Get Timing Parameters
[ opt, oc, tvec, ~, ~, ~ ] = processTimingFixedStep( oc, pulse, opt );
opt = generateRotFrameFrequencyVector( opt );

%% Handle piecewise polynomial order
% RF
if ( ~isfield( pulse, 'orderPP_RF' ) ) && ( ~isfield( pulse, 'orderPP' ) )
    orderPP_RF = uint32( 5 );
    warning( "Set piecewise polynomial order for RF to:\t%i", opt.orderPP_RF );
else
    if ( isfield( pulse, 'orderPP_RF' ) )
        orderPP_RF = uint32( pulse.orderPP_RF );
    elseif ( isfield( pulse, 'orderPP' ) )
        orderPP_RF = uint32( pulse.orderPP );
    end
end
numPPShape_RF = orderPP_RF + 1;

% Grad
if ( ~isfield( pulse, 'orderPP_grad' ) ) && ( ~isfield( pulse, 'orderPP' ) )
    orderPP_grad = uint32( 5 );
    warning( "Set piecewise polynomial order for Grad to:\t%i", opt.orderPP_grad );
else
    if ( isfield( pulse, 'orderPP_grad' ) )
        orderPP_grad = uint32( pulse.orderPP_grad );
    elseif ( isfield( pulse, 'orderPP' ) )
        orderPP_grad = uint32( pulse.orderPP );
    end
end
numPPShape_grad = orderPP_grad + 1;

% Shim
if opt.numZCoils > 0
    if ( ~isfield( pulse, 'orderPP_shim' ) ) && ( ~isfield( pulse, 'orderPP' ) )
        orderPP_shim = uint32( 10 );
        warning( "Set piecewise polynomial order for Shim to:\t%i", opt.orderPP_shim );
    else
        if ( isfield( pulse, 'orderPP_shim' ) )
            orderPP_shim = uint32( pulse.orderPP_shim );
        elseif ( isfield( pulse, 'orderPP' ) )
            orderPP_shim = uint32( pulse.orderPP );
        end
    end

    if ~isfield( opt, 'orderPP_shim' ) || ( opt.orderPP_shim == 0 )
        numPPShape_shim = 0;
        opt.numZCoils = 0;
    else
        numPPShape_shim = orderPP_shim + 1;
    end

else
    orderPP_shim = 0;
    numPPShape_shim = 0;
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

%% Get idxs for the different phases
% RF
[ PPIdxs_RF, numPPPeriods_RF ] = generatePPtvecIdxs( PPStartStop_RF, tvec, idxtol );

% Grad
[ PPIdxs_grad, numPPPeriods_grad ] = generatePPtvecIdxs( PPStartStop_grad, tvec, idxtol );

% Shim
if opt.numZCoils > 0
    [ PPIdxs_shim, numPPPeriods_shim ] = generatePPtvecIdxs( PPStartStop_shim, tvec, idxtol );
end

%% Get the variable and shape function indices for the different waveforms
% RF
[ opt.PPVarIdxs_RF, opt.PPSFIdxs_RF, opt.PPVarSFIdxs_RF, opt.PPVarIdxsCell_RF, opt.PPVarIdxsCoilCell_RF, opt.numVarsPerChannel_RF ] =...
    generatePPVarIdxs( numPPPeriods_RF, numPPShape_RF, opt.numXYCoils, PPStartStop_RF, idxtol );

% Grad
[ opt.PPVarIdxs_grad, opt.PPSFIdxs_grad, opt.PPVarSFIdxs_grad, opt.PPVarIdxsCell_grad, opt.PPVarIdxsCoilCell_grad, opt.numVarsPerChannel_grad ] =...
    generatePPVarIdxs( numPPPeriods_grad, numPPShape_grad, 3, PPStartStop_grad, idxtol );

% Shim
if opt.numZCoils
    [ opt.PPVarIdxs_shim, opt.PPSFIdxs_shim, opt.PPVarSFIdxs_shim, opt.PPVarIdxsCell_shim, opt.PPVarIdxsCoilCell_shim, opt.numVarsPerChannel_shim ] =...
        generatePPVarIdxs( numPPPeriods_shim, numPPShape_shim, opt.numZCoils, PPStartStop_shim, idxtol );
end


%% Add Piecewise polynomial information
% RF
opt.numPPPeriods_RF = numPPPeriods_RF;
opt.PPIdxs_RF = PPIdxs_RF;
opt.PPStartStop_RF = PPStartStop_RF;
opt.orderPP_RF = orderPP_RF;
opt.numPPShape_RF = numPPShape_RF;
opt.timing_RF = timing_RF;

% Grad
opt.numPPPeriods_grad = numPPPeriods_grad;
opt.PPIdxs_grad = PPIdxs_grad;
opt.PPStartStop_grad = PPStartStop_grad;
opt.orderPP_grad = orderPP_grad;
opt.numPPShape_grad = numPPShape_grad;
opt.timing_grad = timing_grad;

% Shim
if opt.numZCoils > 0
    opt.numPPPeriods_shim = numPPPeriods_shim;
    opt.PPIdxs_shim = PPIdxs_shim;
    opt.PPStartStop_shim = PPStartStop_shim;
    opt.orderPP_shim = orderPP_shim;
    opt.numPPShape_shim = numPPShape_shim;
    opt.timing_shim = timing_shim;
else
    opt.numVarsPerChannel_shim = 0;
end

%% Initialize piecewise polynomial samples

% RF
[ opt.varsToTimepoints_RF, opt.varsToChebByPeriods_RF, ~,...
    ~, opt.shapeFnChebCoeffs_RF ] =...
    generatePPChebEvalMatrices(...
    numPPShape_RF, numPPPeriods_RF, opt.PPVarIdxs_RF, opt.PPVarSFIdxs_RF, opt.PPSFIdxs_RF,...
    PPIdxs_RF, PPStartStop_RF, tvec );

% if matches( opt.structtype, 'opt' )
%     opt.shapeFnValsTimePoints_RF_rep = repmat(...
%         reshape( transpose( opt.varsToTimepoints_RF ), [ 1, opt.numVarsPerChannel_RF, opt.numTimePoints ] ),...
%         [ opt.numPos, 1, 1 ] );
% end

% Grad
[ opt.varsToTimepoints_grad, opt.varsToChebByPeriods_grad, ~,...
    ~, opt.shapeFnChebCoeffs_grad ] =...
    generatePPChebEvalMatrices(...
    numPPShape_grad, numPPPeriods_grad, opt.PPVarIdxs_grad, opt.PPVarSFIdxs_grad, opt.PPSFIdxs_grad,...
    PPIdxs_grad, PPStartStop_grad, tvec );

% if matches( opt.structtype, 'opt' )
%     opt.shapeFnValsTimePoints_grad_rep = repmat(...
%         reshape( transpose( opt.varsToTimepoints_grad ), [ 1, opt.numVarsPerChannel_grad, opt.numTimePoints ] ),...
%         [ opt.numPos, 1, 1 ] );
% end

% Shim
if opt.numZCoils > 0
    [ opt.varsToTimepoints_shim, opt.varsToChebByPeriods_shim, ~,...
        ~, opt.shapeFnChebCoeffs_shim ] =...
        generatePPChebEvalMatrices(...
        numPPShape_shim, numPPPeriods_shim, opt.PPVarIdxs_shim, opt.PPVarSFIdxs_shim, opt.PPSFIdxs_shim,...
        PPIdxs_shim, PPStartStop_shim, tvec );
    
    % if matches( opt.structtype, 'opt' )
    %     opt.shapeFnValsTimePoints_shim_rep = repmat(...
    %         reshape( transpose( opt.varsToTimepoints_shim ), [ 1, opt.numVarsPerChannel_shim, opt.numTimePoints ] ),...
    %         [ opt.numPos, 1, 1 ] );
    % end

end

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal", opt.numVarsPerChannel_RF * opt.numXYCoils;...
        "bimag", opt.numVarsPerChannel_RF * opt.numXYCoils;...
        "grad",  opt.numVarsPerChannel_grad * 3;...
        "shim",  opt.numVarsPerChannel_shim * opt.numZCoils;...
        };

    varArray = checkShimVariables( varArray, opt );

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );
    
    opt.estMaxActiveVarsTimeStep = 2 * opt.numXYCoils * numPPShape_RF +...
        numPPShape_grad * 3 + numPPShape_shim * opt.numZCoils;

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

            opt.RFMax_constr = RFMax / sqrt(2);

            for ii = 1:length( RFMaxidx )
                varIdx = varIdxs{ RFMaxidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = RFMax / sqrt( 2 );
            end

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
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = gradMax;
            end

            % still have to implement grad-max constraints in chebyshev basis
            magConstraints( cc ) = false;

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

            % still have to implement grad-max constraints in chebyshev basis
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
    pulse.Z0 = 50; % default impedance is 50 Ohms
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

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPRFMax, 2*opt.numXYCoils*opt.numPPPeriods_RF );

        case { 'RF-slew-rate' }
            
            opt.RFSlewRate_constr = pulse.constraints( 'RF-slew-rate' ); 

            periodDiff_RF = diff( opt.PPStartStop_RF, 1, 2 );
            scDiff_RF = 2 ./ periodDiff_RF;
            opt.D_RF = repmat( getChebDerivativeMatrix( numPPShape_RF ), [ 1, 1, opt.numPPPeriods_RF ] )...
                .* reshape( scDiff_RF, [ 1, 1, opt.numPPPeriods_RF ] );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPRFSlewRate, 2*opt.numXYCoils*opt.numPPPeriods_RF );

        case { 'RF-accel' }

            opt.RFAccel_constr = pulse.constraints( 'RF-accel' );

            periodDiff_RF = diff( opt.PPStartStop_RF, 1, 2 );
            scDiff_RF = 2 ./ periodDiff_RF;
            opt.D_RF = repmat( getChebDerivativeMatrix( numPPShape_RF ), [ 1, 1, opt.numPPPeriods_RF ] )...
                .* reshape( scDiff_RF, [ 1, 1, opt.numPPPeriods_RF ] );
            opt.D2_RF = pagemtimes( opt.D_RF( 1:end-1, 1:end-1, : ), opt.D_RF );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPRFAccel, 2*opt.numXYCoils*opt.numPPPeriods_RF );

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
            opt.adjDFTmtxVarsByTP = opt.adjDFTmtx * opt.varsToTimepoints_RF;

            opt.RFBandwidth_constr = ( RFMax * 1.0e-5 ) * length( opt.adjfvec ) * opt.adjdf;

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPRFBandwidth, opt.numXYCoils );


        case { 'grad-max' }
            
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPGradMax, 3*opt.numPPPeriods_grad );

        case { 'grad-slew-rate' }
            
            opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' ); 
            
            periodDiff_grad = diff( opt.PPStartStop_grad, 1, 2 );
            scDiff_grad = 2 ./ periodDiff_grad;
            opt.D_grad = repmat( getChebDerivativeMatrix( numPPShape_grad ), [ 1, 1, opt.numPPPeriods_grad ] )...
                .* reshape( scDiff_grad, [ 1, 1, opt.numPPPeriods_grad ] );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPGradSlewRate, 3*opt.numPPPeriods_grad );

        case { 'grad-accel' }

            opt.gradAccel_constr = pulse.constraints( 'grad-accel' ); 

            periodDiff_grad = diff( opt.PPStartStop_grad, 1, 2 );
            scDiff_grad = 2 ./ periodDiff_grad;
            opt.D_grad = repmat( getChebDerivativeMatrix( numPPShape_grad ), [ 1, 1, opt.numPPPeriods_grad ] )...
                .* reshape( scDiff_grad, [ 1, 1, opt.numPPPeriods_grad ] );
            opt.D2_grad = pagemtimes( opt.D_grad( 1:end-1, 1:end-1, : ), opt.D_grad );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPGradAccel, 3*opt.numPPPeriods_grad );

        case { 'shim-max' }

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPShimMax, opt.numZCoils*opt.numPPPeriods_shim );

        case { 'shim-slew-rate' }
        
            opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' ); 
            
            periodDiff_shim = diff( opt.PPStartStop_shim, 1, 2 );
            scDiff_shim = 2 ./ periodDiff_shim;
            opt.D_shim = repmat( getChebDerivativeMatrix( numPPShape_shim ), [ 1, 1, opt.numPPPeriods_shim ] )...
                .* reshape( scDiff_shim, [ 1, 1, opt.numPPPeriods_shim ] );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPShimSlewRate, opt.numZCoils*opt.numPPPeriods_shim );

        case { 'shim-accel' }
        
            opt.shimAccel_constr = pulse.constraints( 'shim-accel' ); 
            
            periodDiff_shim = diff( opt.PPStartStop_shim, 1, 2 );
            scDiff_shim = 2 ./ periodDiff_shim;
            opt.D_shim = repmat( getChebDerivativeMatrix( numPPShape_shim ), [ 1, 1, opt.numPPPeriods_shim ] )...
                .* reshape( scDiff_shim, [ 1, 1, opt.numPPPeriods_shim ] );
            opt.D2_shim = pagemtimes( opt.D_shim( 1:end-1, 1:end-1, : ), opt.D_shim );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPShimAccel, opt.numZCoils*opt.numPPPeriods_shim );

        case { 'shim-total' }
        
            opt.shimTotal_constr = pulse.constraints( 'shim-total' ); 
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPShimTotal, numPPPeriods_shim );

        case { 'total-RF-power' }

            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPMaxRFPower, 1 );

        case { 'peak-local-SAR' }

            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPPeakLocalSAR, 1 );

        case { 'peak-global-SAR' }

            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPPeakGlobalSAR, 1 );

        case { 'average-local-SAR' }

            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintPPAvgGlobalSAR, 1 );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generateWaveforms = @generatePPWaveform_base;
opt.generatePlotWaveforms = @generatePPPlotWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientPP_base;

if matches( opt.structtype, 'opt' )
    if opt.useGPU
        opt.gpuArrayAdjointFunction = @gpuArrayAdjointPP_base;
        opt = prepareGPUArrays( opt );

        % opt.gpu_shapeFnValsTimePoints_RF_rep = gpuArray( opt.shapeFnValsTimePoints_RF_rep );
        % opt.gpu_shapeFnValsTimePoints_grad_rep = gpuArray( opt.shapeFnValsTimePoints_grad_rep );
        % if opt.numZCoils > 0
        %     opt.gpu_shapeFnValsTimePoints_shim_rep = gpuArray( opt.shapeFnValsTimePoints_shim_rep );
        % end
    end
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'length' ) = {opt.pulseLength};
pulseinfo( 'numTimePoints' ) = {opt.numTimePoints};
pulseinfo( 'dt' ) = {opt.dt};
pulseinfo( 'dutyCycle' ) = {opt.dutyCycle};
pulseinfo( 'orderPP_RF' ) = {opt.orderPP_RF};
pulseinfo( 'timing_RF' ) = {opt.timing_RF};
pulseinfo( 'orderPP_grad' ) = {opt.orderPP_grad};
pulseinfo( 'timing_grad' ) = {opt.timing_grad};
if opt.numZCoils > 0
    pulseinfo( 'orderPP_shim' ) = {opt.orderPP_shim};
    pulseinfo( 'timing_shim' ) = {opt.timing_shim};
end

pulse.pulseinfo = pulseinfo;

end