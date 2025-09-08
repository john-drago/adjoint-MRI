%% File Initialization
restoredefaultpath;
clear;
close all;
home;

%% File Description
% Base script to optimize selective universal pulses using the adjoint
% method
%
% Specific pulse type: spokes initialization

%% Add Note To Describe Optimization
note = "Universal selective pulse design - spokes initialization";

%% Add paths to directories needed for optimization
currFile = strcat( mfilename( 'fullpath' ), ".m" );
currDir = fileparts( currFile );

cld = strsplit( currDir, 'scripts' );
cld = cld{ 1 };
addpath( genpath( fullfile( cld, 'util' ) ) );
addpath( genpath( fullfile( cld, 'adjoint' ) ) );
tld = normalizePath( fullfile( cld, '..' ) );

dt = datetime('now');
dtIden = dt;
dtIden.Format = 'yyMMddHHmmss';
timestr = strcat('__', char(dtIden) );

fprintf( "\n---------------------------------------------------------\n" );
fprintf( "Script start:\t%s\n", string(dt) );
fprintf( "Script Iden:\t%s", string(dtIden) );
fprintf( "\n---------------------------------------------------------\n" );

%% Set paths
homeDir = getenv( "HOME" );
% matlabDir = fullfile( homeDir, "MATLAB" );
matlabDir = fullfile( homeDir, "matlab" );

% IPOPT
ipoptLibPath = fullfile( matlabDir, "mexIPOPT", "toolbox", "lib" );
ipoptBinPath = fullfile( matlabDir, "Ipopt", "compiled", "bin" );
addpath( genpath( ipoptLibPath ) );
addpath( genpath( ipoptBinPath ) );

%% File control
savePathTop = fullfile( tld, 'data', 'opt', ...
    '2025-05-01-slice-selective-spokes-init',...
    strcat( 'opt', timestr ) );

numWorkers = 1;
binVis = true;
useGPU = true;
saveResult = false;
useParallel = false;
saveFileRecord = true;
ensureFeasibleStart = true;
trackConvergence = true;
trackDecisionVariables = true;

val_di = 5e-3;
withinSlice_di = 0.002;
outsideSlice_di = 0.015;
intraSlice_di = 0.025;
outsideSliceExt = 1.0;
% outsideSlice_di_start = 0.002;
% outsideSlice_di_end = 0.005;

% val_di = 1e-3;
% withinSlice_di = 0.00025;
% outsideSlice_di = 0.002;
% intraSlice_di = 0.010;
% outsideSliceExt = 1.50;
% % outsideSlice_di_start = 0.002;
% % outsideSlice_di_end = 0.005;

%% Name fields and path
subjIden = "pTx_1x8_cylindrical_decoupled";
% subjIden = "BC2_back_feeding";
subjPath = fullfile( tld, 'data', 'fields', '2025-01-01-simulated-fields', subjIden );

fieldsPath = fullfile( strcat( subjPath, ".mat" ) );
VOPpath = fullfile( tld, 'data', 'fields', '2025-01-01-simulated-fields', 'VOPs',...
    strcat( subjIden, '.mat' ) );

% shimPath = fullfile( tld, 'data', 'fields', '2025-01-03-sim-32ch-shim-coil',...
%     '32ch_SensFields.mat');

shimPath = [];

%% Initial Process of fields
UP = load( fieldsPath );
% BS = load( shimPath );

%% Get universal test and train data sets
% numTrain = 10;
% numTest = 3;
% randSplit = randperm( length( UP.subjIden ) );
% 
% universalTrain = UP.subjIden( randSplit( 1:numTrain ) );
% universalTest = UP.subjIden( randSplit( (numTrain+1):(numTrain+numTest) ) );
% 
% idenLogical = contains( UP.subjIden, [ universalTrain; universalTest ], 'ignorecase', true );

idenLogical = true;

%% Fill-in fields struct
[ ~, ~, K ] = ndgrid( 1:length(UP.x), 1:length(UP.y), 1:length(UP.z) );

fields = struct;
fields.x = UP.x;
fields.y = UP.y;
fields.z = UP.z;
fields.X = UP.X;
fields.Y = UP.Y;
fields.Z = UP.Z;
% fields.b1p = UP.b1p( :, :, :, idenLogical ); % for birdcage coil
fields.b1p = UP.b1p( :, :, :, :, idenLogical ); % for pTx coil
fields.db0 = UP.db0shim( :, :, :, idenLogical );
fields.opt_roi = UP.roi_brain( :, :, :, idenLogical );
fields.val_roi = UP.roi_body( :, :, :, idenLogical );

% fields.zi = UP.zi( idenLogical );
% fields.subjIden = string( UP.subjIden( idenLogical ) );
% fields.subjPath = string( UP.subjPath( idenLogical ) );

fields.zi = unique( K( UP.roi_brain ) );
fields.subjIden = subjIden;
fields.subjPath = subjPath;

% fields.bz = zeros( [ BS.coilNum, size( fields.X ) ] );
% for cc = 1:BS.coilNum
%     fields.bz( cc, :, :, : ) = Interp3D(...
%         BS.X, BS.Y, BS.Z, squeeze(BS.bz( cc, :, :, : )),...
%         UP.X, UP.Y, UP.Z );
% end
% clear BS;

fields.bz = 0;

% Save path names
fields.subjectMaps = fieldsPath;
fields.zCoilPath = shimPath;

%% Load VOPs
fields.VOPpath = VOPpath;
VOPst = load( VOPpath );
fields.VOPs = VOPst.VOPs;

% fields.QGlobal = VOPst.QPartialBody;
fields.QGlobal = VOPst.QGlobalMat;

%% Define Pulse Structure
pulse = struct;

constraints = {...
    "total-RF-power", 100;... % W
    "max-RF-power", 24;... % W
    % "RF-slew-rate", 1e7;... % V/s
    % "RF-accel", 1.5e12;... % V/(sec^2)
    "RF-max", 400;... % 2.190253975675163e+02;... % V
    % "peak-local-SAR", 1e4;... % W/kg
    % "peak-global-SAR", 1e4;... % W/kg
    % "RF-bandwidth", 12.5e3;...
    "average-local-SAR", 20;... % W/kg
    "average-global-SAR", 3.2;... % W/kg
    "grad-max", 50e-3;... % T/m
    "grad-slew-rate", 200;... % T/m/sec
    % "grad-accel", 1e8;... % T/m/(sec^2)
    % "shim-max", 30;... % Amp-turns
    % "shim-slew-rate", 3e6;... % Amp-turns/sec
    % "shim-accel", 1e12;... % Amp-turns/(sec^2)
    % "shim-total", 150;... % Amp-turns
    };

if ~isempty( constraints )
    constraintName = string( constraints( :, 1 ) );
    constraintValue = cell2mat( constraints( :, 2 ) );
else
    constraintName = [];
    constraintValue = [];
end

pulse.constraints = dictionary( constraintName, constraintValue );
pulse.excitationType = "slice-selective";
pulse.terminalCostFunction = "FA-with-mean-phase-across-slice";
pulse.runningCostFunction = [];

% specify tailored/universal identity if needed
% pulse.optPopulation = "universal"; % e.g., universal or tailored
% pulse.universalTrain = sort( fields.universalTrain );
% pulse.universalTest = sort( fields.universalTest );

pulse.optPopulation = "tailored"; % e.g., universal or tailored

pulse.dutyCycle = 0.20; % duty cycle in percent
pulse.targFAVal = 10; % degrees
pulse.flipAngleRangePlot = [ 0 20 ];
pulse.dt = 10e-6;
% spokes specific parameters
pulse.numSpokes = 1;
pulse.centralSpokeTBW = 8;
pulse.centralSpokeLength = 1050 * 1e-6;
pulse.nonCentralSpokeTBW = 0;

% pulse.dutyCycle = 0.05; % duty cycle in percent
% pulse.targFAVal = 90; % degrees
% pulse.flipAngleRangePlot = [ 0 120 ];
% pulse.dt = 20e-6;
% % spokes specific parameters
% pulse.numSpokes = 3;
% pulse.centralSpokeTBW = 8;
% pulse.centralSpokeLength = 1240 * 1e-6;
% pulse.nonCentralSpokeTBW = 2;

pulse.targSt = struct;
pulse.targSt.M0Vals = [ 0; 0; 1 ];

pulse.sliceThickness = 5e-3;
pulse.sliceLocation = -6.0 * 1e-2;
pulse.sliceDirection = [ 0; 0; 1 ];

pulse.gyro = 267.5e6;

RFMaxConstr = 0.95 * pulse.constraints( 'RF-max' );
gradSlewConstr = 0.95 * pulse.constraints( 'grad-slew-rate' );
gradMaxConstr = 0.95 * pulse.constraints( 'grad-max' );

roundTime = pulse.dt;
[ spokesTiming ] = determineSpokesDuration( pulse,...
    gradMaxConstr, gradSlewConstr, roundTime );
spokes = determineSpokesTimings( spokesTiming );

pulse.timing_dwxy = spokes.tdom_spokes;
pulse.vals_dwxy = ( -1 * 2*pi ) * spokes.RF_freq;

pulse.centralSpokeTBW = spokesTiming.centralSpokeTBW;
pulse.centralSpokeLength = spokesTiming.centralSpokeLength;
pulse.nonCentralSpokeTBW = spokesTiming.nonCentralSpokeTBW;
% pulse.nonCentralSpokeLength = spokesTiming.nonCentralSpokeLength;
pulse.numSlices = spokesTiming.numSlices;

pulse.constantRotatingFrame = false;
pulse.convertMBackToLarmor = false;
pulse.convertMtargAwayLarmor = false;

pulse.length = spokesTiming.pulseLength;

% % PWC
% pulse.name = "PWC"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
% pulse.type = "base"; % e.g., base or extended for kTP

% % PP
% pulse.name = "PP"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
% pulse.type = "base"; % e.g., base or extended for kTP
% pulse = getPPSpokesTiming( pulse, spokesTiming );
% pulse.orderPP_RF = 20;
% pulse.orderPP_grad = 20;
% pulse.orderPP_shim = 0;
% % pulse.orderPP_RF = 13;
% % pulse.orderPP_grad = 10;
% % pulse.orderPP_shim = 0;

% Cheb
pulse.name = "cheb"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
pulse.type = "base"; % e.g., base or extended for kTP
pulse.orderCheb_RF = 35;
pulse.orderCheb_grad = 50;
pulse.orderCheb_shim = 0;
% pulse.orderCheb_RF = 50;
% pulse.orderCheb_grad = 50;
% pulse.orderCheb_shim = 0;

% % Fourier
% pulse.name = "fourier";
% pulse.type = "base";
% pulse.orderFourier_RF = 17;
% pulse.orderFourier_grad = 25;
% pulse.orderFourier_shim = 0;
% % pulse.orderFourier_RF = 30;
% % pulse.orderFourier_grad = 25;
% % pulse.orderFourier_shim = 0;

% % Fourier-Cheb
% pulse.name = "fourier";
% pulse.type = "cheb";
% pulse.orderFourier_RF = 17;
% pulse.orderCheb_grad = 50;
% pulse.orderCheb_shim = 0;
% % pulse.orderFourier_RF = 30;
% % pulse.orderCheb_grad = 35;
% % pulse.orderCheb_shim = 0;

%% Optimization Control (oc)
oc = struct; % define opt control struct
oc.note = note;
oc.binVis = binVis; % whether or not to make images visible
oc.saveResult = saveResult; % whether or not to save optimization results
oc.saveDir = savePathTop; % location or where to save optimization results

oc.opt_di = struct;

oc.opt_di.withinSlice_di = withinSlice_di;
oc.opt_di.outsideSlice_di = outsideSlice_di;
oc.opt_di.intraSlice_di = intraSlice_di;
oc.opt_di.outsideSliceExt = outsideSliceExt;
% oc.opt_di.outsideSlice_di_start = outsideSlice_di_start;
% oc.opt_di.outsideSlice_di_end = outsideSlice_di_end;

pulse.sliceBounds = determineSliceBounds( pulse.sliceLocation, pulse.sliceThickness );

[ oc.opt_di.z, pulse.sliceBounds ] = getOptimizationDiscretizationLinspaceSliceSelect(...
    withinSlice_di, outsideSlice_di, outsideSliceExt,...
    pulse.sliceLocation, pulse.sliceThickness, pulse.sliceBounds, [ fields.z(1), fields.z(end) ] );

% [ oc.opt_di.z, pulse.sliceBounds ] = getOptimizationDiscretizationLogspaceSliceSelect(...
%     withinSlice_di, outsideSlice_di_start, outsideSlice_di_end,...
%     outsideSliceExt, pulse.sliceLocation, pulse.sliceThickness, pulse.sliceBounds, [ fields.z(1), fields.z(end) ] );

pulse.targSt.x = fields.x(1) : val_di : fields.x(end);
pulse.targSt.y = fields.y(1) : val_di : fields.y(end);
pulse.targSt.z = unique( [...
    fields.z(1) : val_di : fields.z(end),...
    oc.opt_di.z ] );

oc.opt_di.x = fields.x(1) : intraSlice_di : fields.x(end);
oc.opt_di.y = fields.y(1) : intraSlice_di : fields.y(end);

% x_const = [];
% y_const = [];
% z_const = unique( [...
%     pulse.sliceBounds( : );...
%     pulse.sliceBounds( :, 1 ) - withinSlice_di;...
%     pulse.sliceBounds( :, 2 ) + withinSlice_di  ] );
% [ oc.opt_di.x, oc.opt_di.y, oc.opt_di.z ] =  perturbOptimizationPoints(...
%     oc.opt_di.x, oc.opt_di.y, oc.opt_di.z, x_const, y_const, z_const );

oc.opt_dt = pulse.dt; % time in seconds

pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
pulse.targSt = generateMTargSliceSelect( pulse.targSt, pulse );

oc.val_di = zeros( 3, 1 );
oc.val_di( 1:3, 1 ) = val_di;
oc.val_dt = pulse.dt; % time in seconds

maxOptTime = 1 * 60;
maxOptIterations = 300;

% optimization specific parameters

oc.optType = "fmin";
optDisplay = "iter";
oc.fminopt = optimoptions( "fmincon" );
oc.fminopt.Display = optDisplay;
oc.fminopt.UseParallel = useParallel;
oc.fminopt.MaxFunctionEvaluations = inf;
oc.fminopt.SpecifyObjectiveGradient = true;
oc.fminopt.SpecifyConstraintGradient = true;
oc.fminopt.CheckGradients = false;
oc.fminopt.StepTolerance = 1e-10;
oc.fminopt.FiniteDifferenceStepSize = 1e-7;
oc.fminopt.ConstraintTolerance = 1e-8;
oc.fminopt.FunctionTolerance = 1e-6; % doesn't matter for SQP

oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );
oc.fminopt.MaxIterations = maxOptIterations;

% oc.fminopt.Algorithm = 'interior-point';
% % oc.fminopt.BarrierParamUpdate = 'predictor-corrector';
% oc.fminopt.StepTolerance = 1e-10;
% oc.fminopt.HessianApproximation = 'bfgs';
% oc.fminopt.InitBarrierParam = 1e-5;
% % oc.fminopt.EnableFeasibilityMode = true;
% oc.fminopt.SubproblemAlgorithm = 'cg';
% oc.fminopt.TolProjCG = 1e-8;
% oc.fminopt.TolProjCGAbs = 1e-8;

oc.fminopt.Algorithm = 'active-set';
oc.fminopt.RelLineSrchBnd = 2.50e-4;
oc.fminopt.RelLineSrchBndDuration = inf;
% oc.fminopt.TolConSQP = 1e-8;
% oc.fminopt.MaxSQPIter = 5e3;

% oc.optType = 'ipopt';
% oc.ipopt = struct;
% oc.ipopt.options = struct;
% oc.ipopt.options.ipopt = struct;
% oc.ipopt.options.ipopt.print_user_options = 'yes';
% % oc.ipopt.options.ipopt.print_options_documentation = 'yes';
% % oc.ipopt.options.ipopt.print_advanced_options = 'yes';
% oc.ipopt.options.ipopt.print_level = 5;
% oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
% oc.ipopt.options.ipopt.mu_min = 1e-11;
% oc.ipopt.options.ipopt.mu_max = 1e-5;
% oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';
% % oc.ipopt.options.ipopt.adaptive_mu_globalization = 'obj-constr-filter';
% % oc.ipopt.options.ipopt.mu_strategy = 'monotone';
% % oc.ipopt.options.ipopt.mu_init = 1e-6;
% % oc.ipopt.options.ipopt.mu_max_fact = 1000;
% % oc.ipopt.options.ipopt.mu_linear_decrease_factor = 0.2;
% % oc.ipopt.options.ipopt.mu_superlinear_decrease_power = 1.05;
% oc.ipopt.options.ipopt.tol = 1e-6; % desired convergence tolerance
% oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
% oc.ipopt.options.ipopt.dual_inf_tol = 1e-6;
% oc.ipopt.options.ipopt.compl_inf_tol = 1e-6;
% bound_push_frac = 1e-6;
% oc.ipopt.options.ipopt.bound_push = bound_push_frac;
% oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
% slack_bound_push_frac = 1e-2;
% oc.ipopt.options.ipopt.slack_bound_push = slack_bound_push_frac;
% oc.ipopt.options.ipopt.slack_bound_frac = slack_bound_push_frac;
% oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;
% oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
% oc.ipopt.options.ipopt.alpha_min_frac = 0.005;
% oc.ipopt.options.ipopt.max_soc = 10;
% oc.ipopt.options.ipopt.recalc_y = 'yes';
% oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-2;
% oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
% oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
% oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
% oc.ipopt.options.ipopt.linear_solver = 'mumps';
% 
% oc.ipopt.options.ipopt.max_iter = maxOptIterations;
% oc.ipopt.options.ipopt.max_wall_time = maxOptTime;

oc.saveFileRecord = saveFileRecord;
oc.currFile = currFile;
oc.useGPU = useGPU;
oc.trackConvergence = trackConvergence;
oc.trackDecisionVariables = trackDecisionVariables;
oc.ensureFeasibleStart = ensureFeasibleStart;

% Parallel Test
oc = determineParallelWorkers( numWorkers, oc );

%% Initialize opt and val structs
opt = struct;
opt.structtype = 'opt';

timeResDigits = 6; % round to microseconds
opt.di = oc.opt_di;
opt.dt = round( oc.opt_dt, timeResDigits);

% Process fields struct 
% Get data about the fields such as number of coils and number of subjects
[ oc, pulse, fields ] = processFieldsData( oc, pulse, fields );

% Now add fields to the opt struct
[ oc, pulse, fields, opt ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, opt );

% Process Pulse struct
[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );

%% Initialize spokes
p = 8;
spokesFMinIterations = 300;
% spokesFMinIterations = 50;
[ RF0init, K0init, ~, spokes ] = findInitialSpokes( opt, pulse,...
    gradMaxConstr, gradSlewConstr, p, roundTime, spokesFMinIterations );

[ RF0, G0, K0 ] = optimizeSpokesSTA(...
    RF0init, K0init, RFMaxConstr, gradSlewConstr, opt, spokes, spokesFMinIterations );

%% Generate spokes waveform
% RF0init_vec = reshape( RF0init, [ opt.numXYCoils*pulse.numSpokes, 1 ] );
% G0init_vec = reshape( G0init, [ 3*(pulse.numSpokes-1), 1 ] );

% RF0_vec = reshape( RF0, [ opt.numXYCoils*pulse.numSpokes, 1 ] );
% G0_vec = reshape( G0, [ 3*(pulse.numSpokes-1), 1 ] );

opt.RF0 = RF0;
opt.G0 = G0;
opt.K0 = K0;

wvSpokes = makeSpokesWaveformSample( RF0, G0, spokes );

%% Project into PWC basis
% pwccoeffs = getPWCOptCoeffs( wvSpokes, opt.tvec );
% 
% breal_idx_rshp = reshape( opt.breal_idx, [ opt.numTimePoints, opt.numXYCoils ] );
% bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numTimePoints, opt.numXYCoils ] );
% grad_idx_rshp = reshape( opt.grad_idx, [ opt.numTimePoints, 3 ] );
% 
% opt.p0 = zeros( opt.numVars, 1 );
% opt.p0( breal_idx_rshp ) = pwccoeffs.breal;
% opt.p0( bimag_idx_rshp ) = pwccoeffs.bimag;
% opt.p0( grad_idx_rshp ) = pwccoeffs.grad;
% 
% opt.pSc0 = opt.p0 ./ opt.scVec;
% 
% oc.pSc0 = opt.pSc0;
% oc.p0 = opt.p0;

%% Project into PP basis
% ppcoeffs = getPPOptCoeffs( wvSpokes, opt );
% 
% p0 = zeros( opt.numVars, 1 );
% p0( opt.breal_idx ) = reshape( ppcoeffs.breal_coeffs, [ opt.numVarsPerChannel_RF * opt.numXYCoils, 1 ] );
% p0( opt.bimag_idx ) = reshape( ppcoeffs.bimag_coeffs, [ opt.numVarsPerChannel_RF * opt.numXYCoils, 1 ] );
% p0( opt.grad_idx ) = reshape( ppcoeffs.grad_coeffs, [ opt.numVarsPerChannel_grad * 3, 1 ] );
% 
% if opt.numZCoils > 0
%     p0( opt.shim_idx ) = reshape( ppcoeffs.shim_coeffs, [ opt.numVarsPerChannel_shim * opt.numZCoils, 1 ] );
% end
% 
% pSc0 = p0 ./ opt.scVec;
% 
% opt.p0 = p0;
% opt.pSc0 = pSc0;
% 
% oc.pSc0 = opt.pSc0;
% oc.p0 = opt.p0;

%% Project into Cheb basis
chebTol = 1e-6;
chebcoeff = getChebOptCoeffs(...
    wvSpokes, chebTol );

pulse.numCheb_RF = min( size( chebcoeff.breal_coeffs, 1 ), opt.orderCheb_RF + 1 );
pulse.numCheb_grad = min( size( chebcoeff.grad_coeffs, 1 ), opt.orderCheb_grad + 1 );
pulse.numCheb_shim = min( size( chebcoeff.shim_coeffs, 1 ), opt.numCheb_shim );

pulse.orderCheb_RF = opt.numCheb_RF - 1;
pulse.orderCheb_grad = opt.numCheb_grad - 1;
pulse.orderCheb_shim = max( opt.numCheb_shim - 1, 0 );

breal_idx_rshp = reshape( opt.breal_idx, [ opt.numCheb_RF, opt.numXYCoils ] );
bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numCheb_RF, opt.numXYCoils ] );
grad_idx_rshp = reshape( opt.grad_idx, [ opt.numCheb_grad, 3 ] );
if isfield( opt, 'numCheb_shim' ) && ( opt.numCheb_shim > 0 )
    shim_idx_rshp = reshape( opt.shim_idx, [ opt.numCheb_shim, 3 ] );
end

p0 = zeros( opt.numVars, 1 );
p0( breal_idx_rshp ) = chebcoeff.breal_coeffs( 1:pulse.numCheb_RF, : );
p0( bimag_idx_rshp ) = chebcoeff.bimag_coeffs( 1:pulse.numCheb_RF, : );
p0( grad_idx_rshp ) = chebcoeff.grad_coeffs( 1:pulse.numCheb_grad, : );
if isfield( opt, 'numCheb_shim' ) && ( opt.numCheb_shim > 0 )
    p0( shim_idx_rshp ) = chebcoeff.shim_coeffs( 1:pulse.numCheb_shim, : );
end

pSc0 = p0 ./ opt.scVec;

opt.p0 = p0;
opt.pSc0 = pSc0;

oc.pSc0 = opt.pSc0;
oc.p0 = opt.p0;

%% Project into Fourier basis
% fouriercoeff = getFourierOptCoeffs(...
%     wvSpokes, opt.orderFourier_RF, opt.orderFourier_grad, opt.orderFourier_shim );
% 
% breal_idx_rshp = reshape( opt.breal_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
% bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
% grad_idx_rshp = reshape( opt.grad_idx, [ opt.numFourier_grad, 3 ] );
% if isfield( opt, 'numFourier_shim' ) && ( opt.numFourier_shim > 0 )
%     shim_idx_rshp = reshape( opt.shim_idx, [ opt.numFourier_shim, 3 ] );
% end
% 
% p0 = zeros( opt.numVars, 1 );
% p0( breal_idx_rshp ) = fouriercoeff.breal_coeffs;
% p0( bimag_idx_rshp ) = fouriercoeff.bimag_coeffs;
% p0( grad_idx_rshp ) = fouriercoeff.grad_coeffs;
% if isfield( opt, 'numFourier_shim' ) && ( opt.numFourier_shim > 0 )
%     p0( shim_idx_rshp ) = fouriercoeff.shim_coeffs( 1:pulse.numFourier_shim, : );
% end
% 
% pSc0 = p0 ./ opt.scVec;
% 
% opt.p0 = p0;
% opt.pSc0 = pSc0;
% 
% oc.pSc0 = opt.pSc0;
% oc.p0 = opt.p0;

%% Project into Fourier-Cheb basis
% fouriercoeff = getFourierOptCoeffs(...
%     wvSpokes, opt.orderFourier_RF, [], [] );
% 
% breal_idx_rshp = reshape( opt.breal_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
% bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numFourier_RF, opt.numXYCoils ] );
% 
% chebTol = 1e-6;
% chebcoeff = getChebOptCoeffs(...
%     wvSpokes, chebTol );
% 
% pulse.numCheb_grad = min( size( chebcoeff.grad_coeffs, 1 ), opt.orderCheb_grad + 1 );
% pulse.numCheb_shim = min( size( chebcoeff.shim_coeffs, 1 ), opt.numCheb_shim );
% 
% pulse.orderCheb_grad = opt.numCheb_grad - 1;
% pulse.orderCheb_shim = max( opt.numCheb_shim - 1, 0 );
% 
% grad_idx_rshp = reshape( opt.grad_idx, [ opt.numCheb_grad, 3 ] );
% if isfield( opt, 'numCheb_shim' ) && ( opt.numCheb_shim > 0 )
%     shim_idx_rshp = reshape( opt.shim_idx, [ opt.numCheb_shim, 3 ] );
% end
% 
% p0 = zeros( opt.numVars, 1 );
% p0( breal_idx_rshp ) = fouriercoeff.breal_coeffs;
% p0( bimag_idx_rshp ) = fouriercoeff.bimag_coeffs;
% p0( grad_idx_rshp ) = chebcoeff.grad_coeffs( 1:pulse.numCheb_grad, : );
% if isfield( opt, 'numCheb_shim' ) && ( opt.numCheb_shim > 0 )
%     p0( shim_idx_rshp ) = chebcoeff.shim_coeffs( 1:pulse.numCheb_shim, : );
% end
% 
% pSc0 = p0 ./ opt.scVec;
% 
% opt.p0 = p0;
% opt.pSc0 = pSc0;
% 
% oc.pSc0 = opt.pSc0;
% oc.p0 = opt.p0;

%% Prepare for optimization
oc.ensureFeasibleStartRelLineSrchBnd = 5e-5;
oc.ensureFeasibleStartScaleAtEnd = false;

[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
[ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

%% Run adjoint opt
clear costFnTrackDecisionVariablesWrapper;
clear costFnTrackConvergenceWrapper;
clear ipoptObjectiveBestOptWrapper;

 opt = runAdjointOpt( opt, oc );

%% Post process optimization
[ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
    opt, oc, pulse, fields );
