%% File Initialization
% restoredefaultpath;
clear;
% close all;
home;

%% File Description
% Base script to optimize nonselective universal pulses using the adjoint
% method
%
% Specific pulse type:

%% Add Note To Describe Optimization
note = "Tailored selective pulse design from universal";

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

bottomPrevOptDirPath = fullfile( ...
    "2025-03-03-universal-selective",...
    "2025-03-03-universal-spokes-init-90-pwc",...
    "XXXXXXXXXXX",...
    "loc_-2.0cm",...
    "ipopt" );

savePathTop = fullfile( tld, 'data', 'opt',...
    '2025-03-04-universal-2D-excitation-tailored-test',...
    bottomPrevOptDirPath,...
    strcat( 'opt', timestr ) );

prevOptDirPath = fullfile( tld, "data", "opt",...
    bottomPrevOptDirPath );

testSubject = "AdjDataUser112";

%% Set paths
homeDir = getenv( "HOME" );
% matlabDir = fullfile( homeDir, "MATLAB" );
matlabDir = fullfile( homeDir, "matlab" );

% IPOPT
ipoptLibPath = fullfile( matlabDir, "mexIPOPT", "toolbox", "lib" );
ipoptBinPath = fullfile( matlabDir, "Ipopt", "compiled", "bin" );
addpath( genpath( ipoptLibPath ) );
addpath( genpath( ipoptBinPath ) );

%% Initialize Script
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
fieldsPath = fullfile( tld, 'data', 'fields',...
    '2025-01-02-Siemens-Nova-8ch-pTx-Database-scaled', "UPdatabase.mat" );
VOPpath = fullfile( tld, 'data', 'fields',...
    '2025-01-02-Siemens-Nova-8ch-pTx-SAR', "Nova_8ch_pTx_SAR.mat" );

shimPath = [];

% load fields
UP = load( fieldsPath );
VOPst = load( VOPpath );

%% Iterate over the previous optimizations
pattern = '-?\d+(\.\d+)';
sliceLocationInitial = regexp( prevOptDirPath, pattern, 'match'); % Extract number as a string
sliceLocation = str2double( sliceLocationInitial{ end } ) * 1e-2;

locString = sprintf( "loc_%+.1fcm", sliceLocation * 1e2 );

prev = load( fullfile( prevOptDirPath, 'optOutput.mat' ) );

if isfield( prev.oc, 'fminopt' )
    solverName = prev.oc.fminopt.Algorithm;
elseif isfield( prev.oc, 'ipopt' )
    solverName = "ipopt";
else
    error( "Unknown solver name." );
end

fprintf( "\n---------------------------------------------------------\n" );
fprintf( "Time:\t%s\n", string( datetime ) );
fprintf( "In directory:\t%s\n", prevOptDirPath );
fprintf( "Location:\t%+.1fcm\n", sliceLocation * 1e2 );
fprintf( "Test subject:\t%s\n", testSubject );
fprintf( "Using solver:\t%s", solverName );
fprintf( "\n---------------------------------------------------------\n" );

%% Get validation data set
idenLogical = contains( UP.subjIden, testSubject, 'ignorecase', true );

%% Fill-in fields struct
fields = struct;
fields.x = UP.x;
fields.y = UP.y;
fields.z = UP.z;
fields.X = UP.X;
fields.Y = UP.Y;
fields.Z = UP.Z;
fields.b1p = UP.b1p( :, :, :, :, idenLogical );
fields.db0 = UP.db0shim( :, :, :, idenLogical );
fields.opt_roi = UP.roi_brain( :, :, :, idenLogical );
fields.val_roi = UP.roi_body( :, :, :, idenLogical );

fields.zi = UP.zi( idenLogical );
fields.subjIden = string( UP.subjIden( idenLogical ) );
fields.subjPath = string( UP.subjPath( idenLogical ) );

fields.bz = 0;

% Save path names
fields.subjectMaps = fieldsPath;
fields.zCoilPath = shimPath;

%% Load VOPs
fields.VOPpath = VOPpath;
fields.VOPs = VOPst.VOPs;
fields.QGlobal = VOPst.QPartialBody;

%% Define Pulse Structure
pulse = prev.pulse;
pulse.optPopulation = "tailored";
pulse.universalTrain = [];
pulse.universalTest = [];

pulse.dutyCycle = prev.pulse.dutyCycle; % duty cycle in percent
pulse.targFAVal = prev.pulse.targFAVal; % degrees
pulse.flipAngleRangePlot = prev.pulse.flipAngleRangePlot;

% spokes specific parameters
pulse.sliceThickness = prev.pulse.sliceThickness;
pulse.sliceLocation = sliceLocation;
pulse.sliceDirection = [ 0; 0; 1 ];

pulse.numSpokes = prev.pulse.numSpokes;
pulse.centralSpokeTBW = prev.pulse.centralSpokeTBW;
pulse.centralSpokeLength = prev.pulse.centralSpokeLength;
pulse.nonCentralSpokeTBW = prev.pulse.nonCentralSpokeTBW;

pulse.targSt = struct;
pulse.targSt.M0Vals = [ 0; 0; 1 ];

RFMaxConstr = 0.95 * pulse.constraints( 'RF-max' );
gradSlewConstr = 0.95 * pulse.constraints( 'grad-slew-rate' );
gradMaxConstr = 0.95 * pulse.constraints( 'grad-max' );

pulse.dt = prev.pulse.dt;

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

maxOptTime = 20;
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
oc.fminopt.StepTolerance = 1e-8;
oc.fminopt.FiniteDifferenceStepSize = 1e-7;
oc.fminopt.ConstraintTolerance = 1e-6;
oc.fminopt.FunctionTolerance = 1e-6; % doesn't matter for SQP

oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );
oc.fminopt.MaxIterations = maxOptIterations;

oc.fminopt.Algorithm = 'interior-point';
oc.fminopt.BarrierParamUpdate = 'monotone';
oc.fminopt.HessianApproximation = 'bfgs';
oc.fminopt.InitBarrierParam = 1e-5;
oc.fminopt.EnableFeasibilityMode = true;
oc.fminopt.SubproblemAlgorithm = 'cg';
oc.fminopt.StepTolerance = 1e-10;
oc.fminopt.TolProjCG = 1e-8;
oc.fminopt.TolProjCGAbs = 1e-8;

% oc.fminopt.Algorithm = 'active-set';
% oc.fminopt.RelLineSrchBnd = 1.50e-4;
% oc.fminopt.RelLineSrchBndDuration = inf;
% % oc.fminopt.TolConSQP = 1e-6;
% % oc.fminopt.MaxSQPIter = 5e4;

% oc.optType = 'ipopt';
% oc.ipopt = struct;
% oc.ipopt.options = struct;
% oc.ipopt.options.ipopt = struct;
% oc.ipopt.options.ipopt.print_user_options = 'yes';
% % oc.ipopt.options.ipopt.print_options_documentation = 'yes';
% % oc.ipopt.options.ipopt.print_advanced_options = 'yes';
% oc.ipopt.options.ipopt.print_level = 5;
% oc.ipopt.options.ipopt.tol = 1e-4; % desired convergence tolerance
% oc.ipopt.options.ipopt.constr_viol_tol = 1e-4;
% oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
% oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';
% oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
% oc.ipopt.options.ipopt.mu_min = 1e-11;
% oc.ipopt.options.ipopt.mu_max = 1e-5;
% bound_push_frac = 1e-4;
% oc.ipopt.options.ipopt.bound_push = bound_push_frac;
% oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
% oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;
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

%% Prepare for optimization
opt.p0 = prev.opt.pOpt;
opt.pSc0 = prev.opt.pScOpt;

oc.pSc0 = opt.pSc0;
oc.p0 = opt.p0;

[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );

if saveResult && strcmpi( solverName, "ipopt" )
    if ~isfolder( fullfile( saveDir, testSubject ) )
        mkdir( fullfile( saveDir, testSubject ) )
    end
    oc.ipopt.options.ipopt.output_file = char( fullfile( oc.saveDir, testSubject, "IPOPT_OUTPUT.txt" ) );
end
[ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

%% Run adjoint opt
clear costFnTrackDecisionVariablesWrapper;
clear costFnTrackConvergenceWrapper;
clear ipoptObjectiveBestOptWrapper;

opt = runAdjointOpt( opt, oc );

%% Post process
oc.saveDir = fullfile( oc.saveDir, testSubject );
[ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
    opt, oc, pulse, fields );