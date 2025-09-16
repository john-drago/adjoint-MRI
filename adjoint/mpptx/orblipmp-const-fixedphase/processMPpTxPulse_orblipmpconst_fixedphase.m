function [ oc, pulse, opt ] = processMPpTxPulse_orblipmpconst_fixedphase( oc, pulse, opt )
% This function will initialize an "or-blip-mp" mp-ptx pulse and generate
% the needed constraints

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Get Timing Parameters
[ opt, oc, ~, ~, ~, ~ ] = processTimingFixedStep( oc, pulse, opt );
opt = generateRotFrameFrequencyVector( opt );

%% Determine MPpTx Timing
[ oc, pulse, opt ] = processMPpTxTiming( oc, pulse, opt );

%% Add slew rate information to opt struct
if opt.numZCoils > 0
    opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
end
opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal-ORSP", opt.numXYCoils;...
        "bimag-ORSP", opt.numXYCoils;...
        "breal-MPSP", opt.numXYCoils;...
        "bimag-MPSP", opt.numXYCoils;...
        "Gx-Blip", ( 1 );...
        "Gy-Blip", ( 1 );...
        "Gz-Blip", ( 1 );...
        "Gxreal-MPSP", ( 1 );...
        "Gyreal-MPSP", ( 1 );...
        "Gzreal-MPSP", ( 1 );...
        "Gximag-MPSP", ( 1 );...
        "Gyimag-MPSP", ( 1 );...
        "Gzimag-MPSP", ( 1 );...
        "shimmag-MPSP", (opt.numZCoils);...
        "shimph-MPSP", 1;...
        "shimblipscale", 1;...
        };

    varNames = string( varArray( :, 1 ) );

    if opt.ORSP_i == opt.ORSP_f
        ORSPVars = contains(varNames, "ORSP", 'ignorecase', true );
        varArray = varArray( ~ORSPVars, : );
    end

    if isempty(opt.Blip_i) && isempty(opt.Blip_f)
        BlipVars = contains(varNames, "Blip", 'ignorecase', true );
        varArray = varArray( ~BlipVars, : );
    end

    if ( isempty(opt.MPSP_i) && isempty(opt.MPSP_f) ) || ( opt.MPSP_i == opt.MPSP_f )
        MPSPVars = contains(varNames, "MPSP", 'ignorecase', true );
        varArray = varArray( ~MPSPVars, : );
    end

    if opt.numZCoils == 0
        shimVars = contains(varNames, "shim", 'ignorecase', true );
        varArray = varArray( ~shimVars, : );
    end

    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );

    opt.estMaxActiveVarsTimeStep = 2 * opt.numXYCoils + 2 * 3 + 1 + opt.numZCoils;

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
        case { 'grad-max-MPSP' }
            gradMax = pulse.constraints( 'grad-max-MPSP' );

            opt.gradMaxMPSP_constr = gradMax;

            % Then handle MPSP gradient constraints
            gradMaxMPSPidx = find( contains( varNames,...
                [ "Gxreal"; "Gyreal"; "Gzreal"; "Gximag"; "Gyimag"; "Gzimag" ],...
                'ignorecase', true ) );
            for ii = 1:length( gradMaxMPSPidx )
                varIdx = varIdxs{ gradMaxMPSPidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = gradMax / sqrt(2);
            end
        case { 'grad-blip-max' }
            gradBlipMax = pulse.constraints( 'grad-blip-max' );

            opt.gradBlipMax_constr = gradBlipMax;

            % Handle blip constraints
            gradMaxBlipidx = find( contains( varNames,...
                [ "Gx-Blip"; "Gy-Blip"; "Gz-Blip" ],...
                'ignorecase', true ) );
            for ii = 1:length( gradMaxBlipidx )
                varIdx = varIdxs{ gradMaxBlipidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = gradBlipMax;
            end
        case { 'shim-max-MPSP' }
            shimMax = pulse.constraints( 'shim-max-MPSP' );
            opt.shimMaxMPSP_constr = shimMax;
            % Handle MPSP gradient constraints
            shimMaxMPSPidx = find( contains( varNames,...
                "shimmag",...
                'ignorecase', true ) );
            for ii = 1:length( shimMaxMPSPidx )
                varIdx = varIdxs{ shimMaxMPSPidx( ii ) };
                lb( varIdx ) = -1;
                ub( varIdx ) = 1;
                scVec( varIdx ) = shimMax;
            end
        case { 'shim-blip-max' }
            shimBlipMax = pulse.constraints( 'shim-blip-max' );
            opt.shimBlipMax_constr = shimBlipMax;
            % Handle blip constraints
            shimMaxBlipidx = find( contains( varNames,...
                "shimblip",...
                'ignorecase', true ) );
            shimMax = pulse.constraints( 'shim-max-MPSP' );
            % Need ratio of shim Blip to shim Max MPSP
            varIdx = varIdxs{ shimMaxBlipidx }; %#ok
            lb( varIdx ) = -1;
            ub( varIdx ) = 1;
            scVec( varIdx ) = shimBlipMax / shimMax;
        case { 'RF-slew-rate' }
            continue % already took care of this
        case { 'grad-slew-rate' }
            continue % already took care of this
        case { 'grad-max' }
            continue % already took care of this
        case { 'shim-slew-rate' }
            continue
        case { 'shim-max' }
            continue
        otherwise
            % assume constraint is specified linear inequality or convex constraint and pass to
            % next stage. If error, will catch during next stage.
            magConstraints( cc ) = false;
            continue
    end
end

if opt.numZCoils > 0
    % Scale phase with 2pi
    scVec( varIdxs{ contains( varNames, "shimph", 'ignorecase', true ) } ) = 2*pi;
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
        case { 'peak-local-SAR' }

            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulsePeakLocalSAR, 1 );

        case { 'peak-global-SAR' }
            
            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulsePeakGlobalSAR, 1 );

        case { 'average-local-SAR' }
            
            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulseAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulseAvgGlobalSAR, 1 );

        case { 'total-RF-power' }
            
            opt.totalRFPower_constr = pulse.constraints( "total-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulseTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxConstPulseMaxRFPower, 1 );

        case { 'shim-total' }
            
            opt.shimTotal_constr = pulse.constraints( "shim-total" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintMPpTxPulseFixedPhaseShimTotal, 1 );

        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generatePlotWaveforms = @generateMPpTxPlotWaveform_orblipmpconst_fixedphase;
opt.generateWaveforms = @generateMPpTxWaveform_orblipmpconst_fixedphase;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientMPpTx_orblipmpconst_fixedphase;

if opt.useGPU
    opt.gpuArrayAdjointFunction = @gpuArrayAdjointMPpTx_orblipmpconst_fixedphase;
    opt = prepareGPUArrays( opt );
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'fz_mpptx' ) = opt.fz_mpptx;
pulseinfo( 'wz_mpptx' ) = opt.wz_mpptx;
pulseinfo( 'dfxy_mpptx' ) = opt.dfxy_mpptx;
pulseinfo( 'dwxy_mpptx' ) = opt.dwxy_mpptx;

pulseinfo( 'tORSP' ) = pulse.tORSP;
pulseinfo( 'tBlip' ) = pulse.tBlip;
pulseinfo( 'tMPSP' ) = pulse.tMPSP;

pulseinfo( 'minSlewTime' ) = opt.minSlewTime;
pulseinfo( 'minRFSlewTime' ) = opt.minRFSlewTime;

pulse.pulseinfo = pulseinfo;

end
