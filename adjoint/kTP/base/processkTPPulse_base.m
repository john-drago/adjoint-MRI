function [ oc, pulse, opt ] = processkTPPulse_base( oc, pulse, opt )
% This function will initialize an "or-blip-mp" mp-ptx pulse and generate
% the needed constraints

%% Add kTP parameters to opt struct
opt.num_kTP = pulse.num_kTP;

opt.gyro = 267.5e6;

%% Assign GPU information
[ oc, opt ] = assignGPUInformation( oc, opt );

%% Get Initial Timing Parameters
% timeResDigits = 6; % round to microseconds
% pulseLength = round( pulse.length, timeResDigits); 
dt_min = opt.dt; % get spacing of points for optimization
dt_tol = 1e-6;

% We are adding in the idea that dt could be variable for future pulse
% design opportunities.

%% Determine if the kTP parameters are valid
num_kTP = pulse.num_kTP;
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
        pulse.constraints( "shim-max" ) = 50; % Amp-turns
    end
    if ~any( contains( constraintList, "shim-slew-rate", 'ignorecase', true ) )
        pulse.constraints( "shim-slew-rate" ) = 3e6; % Amp-turns/sec
    end
end

%% Determine slew limits
numSigDigitsRound = 6;

% need to determine the slew time for the RF during the RF periods of kTP
if isfield( pulse, "minRFSlewTime" )
    minRFSlewTime = pulse.minRFSlewTime;
else
    minRFSlewTime = 10e-6; % 10 us minimum slew time
end

if ~isfield( pulse, 'RFSlewTime' ) && any( contains( keys( pulse.constraints ), 'RF-slew-rate', 'IgnoreCase', true ) )
    RFSlewTime = max( [ abs( pulse.constraints( 'RF-max' ) / pulse.constraints( 'RF-slew-rate' ) ), minRFSlewTime ] );
    RFSlewTime = round( round( RFSlewTime/dt_min ) * dt_min, numSigDigitsRound, 'significant');
else
    RFSlewTime = pulse.minRFSlewTime;
end

% need to determine the maximum grad value during the blip period of kTP to
% determine if the constraint will be violated
if isfield( pulse, "minGradSlewTime" )
    minGradSlewTime = pulse.minGradSlewTime;
else
    minGradSlewTime = 10e-6; % 10 us minimum slew time
end

if ~isfield( pulse, 'gradSlewTime' )
    gradSlewTime = max( [ abs( pulse.constraints( 'grad-max' ) / pulse.constraints( 'grad-slew-rate' ) ), minGradSlewTime ] );
    gradSlewTime = round( round( gradSlewTime/dt_min ) * dt_min, numSigDigitsRound, 'significant');
else
    gradSlewTime = pulse.gradSlewTime;
end

% need to determine the maximum shim value during the blip period of kTP to
% determine if the constraint will be violated
if isfield( pulse, "minGradSlewTime" )
    minShimSlewTime = pulse.minShimSlewTime;
else
    minShimSlewTime = 10e-6; % 10 us minimum slew time
end

if opt.numZCoils > 0
    if ~isfield( pulse, 'shimSlewTime' )
        shimSlewTime = max( [ abs( pulse.constraints( 'shim-max' ) / pulse.constraints( 'shim-slew-rate' ) ), minShimSlewTime ] );
        shimSlewTime = round( round( shimSlewTime/dt_min ) * dt_min, numSigDigitsRound, 'significant');
    else
        shimSlewTime = pulse.shimSlewTime;
    end
end

%% Determine Number of kT-points and timing

RFLength = max( [ pulse.RFLength, 2*RFSlewTime ] );
if opt.numZCoils > 0
    blipLength = max( [ pulse.blipLength, 2*gradSlewTime, 2*shimSlewTime ] );
else
    blipLength = max( [ pulse.blipLength, 2*gradSlewTime ] );
end

if RFLength > pulse.RFLength
    warning( "Elongated RF slew length to ensure slew rate met." )
end

if blipLength > pulse.blipLength
    warning( "Elongated blip length to ensure slew rate met." )
end

% Check to ensure that the new timings are valid and if not, adjust them
estkTPlength = num_kTP * RFLength + ( num_kTP - 1 ) * blipLength;
if abs( estkTPlength - pulse.length ) / pulse.length > dt_tol
    kTPlengthWarning = "Have to elongate pulse length to ensure slew rates are satisfied.\n" + ...
        "Consider changing RFLength and blipLength timings to accommodate desired pulse length.\n";
    warning( kTPlengthWarning );
    pulse.length = estkTPlength;
    pulseLength = estkTPlength;
else
    pulseLength = pulse.length;
end

% Assign timings
[ RFTiming, RFCenters, blipTiming, blipCenters ] = ...
    getRFBlipTimings( num_kTP, RFLength, blipLength );

%% Redo timing parameters to account for possibility that pulse elongated
tperiods_unsrt = [ RFCenters; blipCenters ];
[ tperiods, tvec_srt_idx ] = sort( tperiods_unsrt );

n = 1;
dim = 2;
dtperiods_unsrt = diff( [ RFTiming; blipTiming ], n, dim  );
dtperiods = dtperiods_unsrt( tvec_srt_idx );

numOptTimePoints = length( tperiods );
numTimePoints = numOptTimePoints + 2 * num_kTP;

time_idx = 1:numTimePoints;

RF_idx = time_idx( 2:4:end );
RF_Slew_i_idx = time_idx( 1:4:end );
RF_Slew_f_idx = time_idx( 3:4:end );
blip_idx = time_idx( 4:4:end );

dtvec = zeros( numTimePoints, 1 ); 
dtvec( RF_idx ) = dtperiods( 1:2:end ) - 2*RFSlewTime;
dtvec( RF_Slew_i_idx ) = RFSlewTime;
dtvec( RF_Slew_f_idx ) = RFSlewTime;
dtvec( blip_idx ) = dtperiods( 2:2:end );

tvec = cumsum( dtvec ); 
tvec = [ 0; tvec( 1:(end-1) ) ] + dtvec/2;

tvec = transpose( tvec );
dtvec = transpose( dtvec );

if strcmpi( opt.structtype, "val" )
    opt.opt_tvec = oc.opt_tvec;
    opt.opt_dtvec = oc.opt_dtvec;
elseif strcmpi( opt.structtype, "opt" )
    opt.opt_tvec = tvec;
    opt.opt_dtvec = dtvec;
    oc.opt_tvec = tvec;
    oc.opt_dtvec = dtvec;
end

%% Add slew rate information to opt struct
if opt.numZCoils > 0
    opt.shimSlewRate_constr = pulse.constraints( 'shim-slew-rate' );
end
opt.gradSlewRate_constr = pulse.constraints( 'grad-slew-rate' );

%% Remove slew rate constraints because they are accounted for implicitly by timing and max values
% pulse.constraints( 'RF-slew-rate' ) = [];
% pulse.constraints( 'grad-slew-rate' ) = [];
% if opt.numZCoils > 0
%     pulse.constraints( 'shim-slew-rate' ) = [];
% end

%% Assign to opt struct
opt.RFTiming = RFTiming;
opt.RFCenters = RFCenters;
opt.RFSlewTime = RFSlewTime;
opt.RFLength = RFLength;

opt.blipTiming = blipTiming;
opt.blipLength = blipLength;
opt.blipCenters = blipCenters;

opt.gradSlewTime = gradSlewTime;

if opt.numZCoils > 0
    opt.shimSlewTime = shimSlewTime;
end

opt.RF_idx = RF_idx;
opt.RF_Slew_i_idx = RF_Slew_i_idx;
opt.RF_Slew_f_idx = RF_Slew_f_idx;
opt.blip_idx = blip_idx;

opt.numTimePoints = numTimePoints;
opt.pulseLength = pulseLength;
opt.tvec = tvec;
opt.dtvec = dtvec;

opt = generateRotFrameFrequencyVector( opt );

opt.get_tvec = @get_tvec;

opt.num_kTP = num_kTP;

%% Determine optimization parameters
if strcmpi( opt.structtype, "opt" )
    varArray = {...
        "breal", num_kTP * opt.numXYCoils;...
        "bimag", num_kTP * opt.numXYCoils;...
        "grad", 3 * ( num_kTP - 1 );...
        "shim", opt.numZCoils * ( num_kTP - 1 );...
        };

    varArray = checkShimVariables( varArray, opt );
    
    [ opt, ~, varNames, ~, ~, ~, varIdxs, numVars ] =...
        processVarOrganization( opt, varArray );
    
    opt.estMaxActiveVarsTimeStep = max( [ 2 * opt.numXYCoils, 3 + opt.numZCoils ] );

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
        case { 'RF-slew-rate' }
            continue; % already addressed above
        case { 'grad-slew-rate' }
            continue; % already addressed above
        case { 'shim-slew-rate' }
            continue; % already addressed above
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
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulseTotalRFPower, 1 );

        case { 'max-RF-power' }

            opt.maxRFPower_constr = pulse.constraints( "max-RF-power" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulseMaxRFPower, 1 );

        case { 'peak-local-SAR' }
            
            opt.peakLocalSAR_constr = pulse.constraints( "peak-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulsePeakLocalSAR, 1 );

        case { 'peak-global-SAR' }
            
            opt.peakGlobalSAR_constr = pulse.constraints( "peak-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulsePeakGlobalSAR, 1 );

        case { 'average-local-SAR' }

            opt.avgLocalSAR_constr = pulse.constraints( "average-local-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulseAvgLocalSAR, opt.numVOPs );

        case { 'average-global-SAR' }

            opt.avgGlobalSAR_constr = pulse.constraints( "average-global-SAR" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulseAvgGlobalSAR, 1 );

        case { 'shim-total' }

            opt.shimTotal_constr = pulse.constraints( "shim-total" );
            [ nlconIneqSt ] = nlconAddList( nlconIneqSt, @constraintkTPPulseShimTotal, num_kTP - 1 );
            
        otherwise
            error( "Unknown constraints type:\t%s", convConstraintList( cc ) )
    end
end

%% Gather constraints and scaling
opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt );

%% Add function to generate waveforms
opt.generatePlotWaveforms = @generatekTPPlotWaveform_base;
opt.generateWaveforms = @generatekTPWaveform_base;

%% Add function to calculate gradp_f array
opt.forwardModelGradientFunction = @forwardModelGradientkTP_base;

if opt.useGPU
    opt.gpuArrayAdjointFunction = @gpuArrayAdjointkTP_base;
    opt = prepareGPUArrays( opt );
end

%% Add post process function
opt.postProcessFunction = @postProcessAdjointShim;

%% Add pulse information for saving
pulseinfo = dictionary;
pulseinfo( 'num_kTP' ) = num_kTP;
pulseinfo( 'RFLength' ) = RFLength;
pulseinfo( 'blipLength' ) = blipLength;
pulseinfo( 'numTimePoints' ) = numTimePoints;

pulseinfo( 'minRFSlewTime' ) = minRFSlewTime;
pulseinfo( 'minGradSlewTime' ) = minGradSlewTime;
pulseinfo( 'minShimSlewTime' ) = minShimSlewTime;

pulse.pulseinfo = pulseinfo;

pulse.RFLength = RFLength;
pulse.blipLength = blipLength;
pulse.minRFSlewTime = minRFSlewTime;
pulse.minGradSlewTime = minGradSlewTime;
pulse.minShimSlewTime = minShimSlewTime;

end

%% Helper Function
% ----------------------------------------------------------------------- %
function opt = get_tvec( ~, opt )

if strcmpi (opt.structtype, 'opt')
    error( "This should not be called for kTP." )
elseif strcmpi( opt.structtype, 'val' )
    
    tol = 1e-9;

    pulseLength = opt.pulseLength;

    tvec = ( opt.dt/2 : opt.dt : (pulseLength - tol ) );
    dtvec = opt.dt * ones( size( tvec ) );

    opt.tvec = tvec;
    opt.dtvec = dtvec;

    opt.numTimePoints = length( dtvec );

end

end
% ----------------------------------------------------------------------- %