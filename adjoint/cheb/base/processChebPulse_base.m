function [ oc, pulse, opt ] = processChebPulse_base( oc, pulse, opt )
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

%% Handle chebyshev polynomials
if isfield( pulse, 'orderCheb_RF' )
    opt.orderCheb_RF = uint32( pulse.orderCheb_RF );
elseif isfield( pulse, 'orderCheb' )
    opt.orderCheb_RF = uint32( pulse.orderCheb );
else
    error( "Did not specify order of RF Cheb polynomial" );
end
opt.numCheb_RF = opt.orderCheb_RF + 1;

if isfield( pulse, 'orderCheb_grad' )
    opt.orderCheb_grad = uint32( pulse.orderCheb_grad );
elseif isfield( pulse, 'orderCheb' )
    opt.orderCheb_grad = uint32( pulse.orderCheb );
else
    error( "Did not specify order of grad Cheb polynomial" );
end
opt.numCheb_grad = opt.orderCheb_grad + 1;

if opt.numZCoils > 0
    if isfield( pulse, 'orderCheb_shim' )
        opt.orderCheb_shim = uint32( pulse.orderCheb_shim );
    elseif isfield( pulse, 'orderCheb' )
        opt.orderCheb_shim = uint32( pulse.orderCheb );
    else
        error( "Did not specify order of shim Cheb polynomial" );
    end
    if ~isfield( opt, 'orderCheb_shim' ) || ( opt.orderCheb_shim == 0 ) 
        opt.numCheb_shim = 0;
        opt.numZCoils = 0;
    else
        opt.numCheb_shim = opt.orderCheb_shim + 1;
    end
    
else
    if isfield( pulse, 'orderCheb_shim' )
        if pulse.orderCheb_shim > 0
            warning( "Removing shim optimization variables, because detected numZCoils = 0." );
        end
    end
    opt.orderCheb_shim = 0;
    opt.numCheb_shim = 0;
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
if ( opt.numZCoils > 0 ) && ( opt.numCheb_shim > 0 )
    if ~any( contains( constraintList, "shim-max", 'ignorecase', true ) )
        pulse.constraints( "shim-max" ) = 30; % Amp/turn
    end
    if ~any( contains( constraintList, "shim-slew-rate", 'ignorecase', true ) )
        pulse.constraints( "shim-slew-rate" ) = 3e6; % Amp/turn/sec
    end
end

%% Add slew rate information to opt struct
opt.gradSlewRate = pulse.constraints( 'grad-slew-rate' );
if ( opt.numZCoils > 0 ) && ( opt.numCheb_shim > 0 )
    opt.shimSlewRate = pulse.constraints( 'shim-slew-rate' );
end

%% Initialize Cheb polynomial samples
if ( opt.numZCoils > 0 ) && ( opt.numCheb_shim > 0 )
    opt.numCheb_max = max( [ opt.numCheb_RF, opt.numCheb_grad, opt.numCheb_shim ] );
else
    opt.numCheb_max = max( [ opt.numCheb_RF, opt.numCheb_grad ] );
end
opt.orderCheb_max = opt.numCheb_max - 1; 

opt.Tn = evalChebClenshaw( transpose( opt.tvec ), eye( opt.numCheb_max ), opt.tdom );

% if matches( opt.structtype, 'opt' )
%     opt.Tn_rep = repmat( reshape( opt.Tn, [ 1, opt.numCheb_max, opt.numTimePoints ] ), [ opt.numPos, 1, 1 ] );
% end

%% Assign Functions to evaluate function and derivative of function
% opt.calcWaveforms = @evalChebClenshaw;

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal", opt.numCheb_RF * opt.numXYCoils;...
        "bimag", opt.numCheb_RF * opt.numXYCoils;...
        "grad",  opt.numCheb_grad * 3;...
        "shim",  opt.numCheb_shim * opt.numZCoils;...
        };

    varArray = checkShimVariables( varArray, opt );

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );
    
    opt.estMaxActiveVarsTimeStep = 2 * opt.numXYCoils * opt.numCheb_RF...
        + opt.numCheb_grad * 3 + opt.numCheb_shim * opt.numZCoils;

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
                % lb( varIdx ) = -1;
                % ub( varIdx ) = 1;
                scVec( varIdx ) = RFMax / sqrt( 2 );
            end

            % Set value at t=0 and t=pulseLength to be zero
            % breal
            [ Aeqbrealif, beqbrealif ] = constraintChebWaveformInitFinal(...
                opt.breal_idx, opt.numCheb_RF, opt );
            AbeqSt = AbAddList( AbeqSt, "breal-if", Aeqbrealif, beqbrealif );
            
            % bimag
            [ Aeqimagif, beqbimagif ] = constraintChebWaveformInitFinal(...
                opt.bimag_idx, opt.numCheb_RF, opt );
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
            % grad
            [ Aeqgradif, beqgradif ] = constraintChebWaveformInitFinal(...
                opt.grad_idx, opt.numCheb_grad, opt );
            AbeqSt = AbAddList( AbeqSt, "grad-if", Aeqgradif, beqgradif );

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
                % lb( varIdx ) = -1;
                % ub( varIdx ) = 1;
                scVec( varIdx ) = shimMax;
            end

            % Set value at t=0 and t=pulseLength to be zero
            % shim
            [ Aeqshimif, beqshimif ] = constraintChebWaveformInitFinal(...
                opt.shim_idx, opt.numCheb_shim, opt );
            AbeqSt = AbAddList( AbeqSt, "shim-if", Aeqshimif, beqshimif );

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

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebRFMax, 2*opt.numXYCoils );

        case { 'RF-slew-rate' }
            
            opt.RFSlewRate_constr = pulse.constraints( 'RF-slew-rate' ); 
            opt.D_RF = getChebDerivativeMatrix( opt.numCheb_RF, opt.tdom );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebRFSlewRate, 2*opt.numXYCoils );

        case { 'RF-accel' }

            opt.RFAccel_constr = pulse.constraints( 'RF-accel' ); 
            opt.D_RF = getChebDerivativeMatrix( opt.numCheb_RF, opt.tdom );
            opt.D2_RF = opt.D_RF( 1:end-1, 1:end-1 ) * opt.D_RF;
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebRFAccel, 2*opt.numXYCoils );
        
        case { 'grad-max' }
            
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebGradMax, 3 );

        case { 'grad-slew-rate' }

            opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );
            opt.D_grad = getChebDerivativeMatrix( opt.numCheb_grad, opt.tdom );
            maxExtrema = 10;
            maxExtrema = min( [ maxExtrema, opt.orderCheb_grad - 1 ] );
            maxExtrema = max( [ 1, maxExtrema ] );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @( pSc, opt ) constraintChebGradSlewRate( pSc, opt, maxExtrema ), 3 * maxExtrema );

            % [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebGradSlewRate, 3 * opt.numTimePoints );

        case { 'grad-accel' }

            opt.gradAccel_constr = pulse.constraints( 'grad-accel' );
            opt.D_grad = getChebDerivativeMatrix( opt.numCheb_grad, opt.tdom );
            opt.D2_grad = opt.D_grad( 1:end-1, 1:end-1 ) * opt.D_grad;
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebGradAccel, 3 );

        case { 'shim-max' }

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebShimMax, opt.numZCoils );

        case { 'shim-slew-rate' }

            opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
            opt.D_shim = getChebDerivativeMatrix( opt.numCheb_shim, opt.tdom );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebShimSlewRate, opt.numZCoils );

        case { 'shim-accel' }

            opt.shimAccel_constr = pulse.constraints( 'shim-accel' );
            opt.D_shim = getChebDerivativeMatrix( opt.numCheb_shim, opt.tdom );
            opt.D2_shim = opt.D_shim( 1:end-1, 1:end-1 ) * opt.D_shim;
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebShimAccel, opt.numZCoils );

        case { 'shim-total' }

            opt.shimTotal_constr = pulse.constraints( 'shim-total' );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebShimTotal, 1 );

        case { 'total-RF-power' }

            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebMaxRFPower, 1 );

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
            opt.adjDFTmtxTn = opt.adjDFTmtx * opt.Tn( :, 1:opt.numCheb_RF );
            
            opt.RFBandwidth_constr = ( RFMax * 1.0e-5 ) * length( opt.adjfvec ) * opt.adjdf;

            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebRFBandwidth, opt.numXYCoils );

        case { 'peak-local-SAR' }

            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebPeakLocalSAR, 1 );

        case { 'peak-global-SAR' }

            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebPeakGlobalSAR, 1 );

        case { 'average-local-SAR' }

            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintChebAvgGlobalSAR, 1 );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generateWaveforms = @generateChebWaveform_base;
opt.generatePlotWaveforms = @generateChebPlotWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientCheb_base;

if matches( opt.structtype, 'opt' )
    if opt.useGPU
        opt.gpuArrayAdjointFunction = @gpuArrayAdjointCheb_base;
        opt = prepareGPUArrays( opt );
        opt.gpu_Tn = gpuArray( opt.Tn );
        % opt.gpu_Tn_rep = gpuArray( opt.Tn_rep );
    end
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'numTimePoints' ) = numTimePoints;
pulseinfo( 'dt' ) = opt.dt;
pulseinfo( 'dutyCycle' ) = opt.dutyCycle;
pulseinfo( 'orderCheb_max' ) = opt.orderCheb_max;
pulseinfo( 'orderCheb_RF' ) = opt.orderCheb_RF;
pulseinfo( 'orderCheb_grad' ) = opt.orderCheb_grad;
if ( opt.numZCoils > 0 ) && ( opt.numCheb_shim > 0 )
    pulseinfo( 'orderCheb_shim' ) = opt.orderCheb_shim;
else
    pulseinfo( 'orderCheb_shim' ) = 0;
end

pulse.pulseinfo = pulseinfo;

end