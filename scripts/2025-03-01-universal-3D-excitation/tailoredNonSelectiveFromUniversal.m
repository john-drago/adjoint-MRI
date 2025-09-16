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
note = "Tailored nonselective pulse design from universal";

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
    "2025-03-01-universal-nonselective",...
    "2025-03-01-universal-kTP-init-90-cheb",...
    "XXXXXXXXXXXX",...
    "ipopt" );

savePathTop = fullfile( tld, 'data', 'opt',...
    '2025-03-01-universal-3D-excitation-tailored-test',...
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

opt_di = 5e-3;
val_di = 2e-3;

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
fprintf( "Using solver:\t%s\n", solverName );
fprintf( "Test subject:\t%s", testSubject );
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

%% Optimization Control (oc)
oc = struct; % define opt control struct
oc.binVis = binVis; % whether or not to make images visible
oc.saveResult = saveResult; % whether or not to save optimization results
oc.saveDir = savePathTop; % location or where to save optimization results
oc.note = note;

oc.opt_di = opt_di * ones( 3, 1 ); % dx, dy, dz for optimization
oc.opt_dt = prev.oc.opt_dt; % time in seconds

oc.val_di = val_di * ones( 3, 1 ); % dx, dy, dz for validation
oc.val_dt = prev.oc.val_dt; % time in seconds

maxOptTime = 180;
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

oc.fminopt.Algorithm = 'active-set';
oc.fminopt.StepTolerance = 1e-6;
oc.fminopt.RelLineSrchBnd = 1e-3;
oc.fminopt.RelLineSrchBndDuration = inf;
oc.fminopt.TolConSQP = 1e-8;

% oc.fminopt.Algorithm = 'interior-point';
% oc.fminopt.StepTolerance = 1e-10;
% % oc.fminopt.BarrierParamUpdate = 'predictor-corrector';
% oc.fminopt.HessianApproximation = 'bfgs';
% oc.fminopt.InitBarrierParam = 1e-5;
% % oc.fminopt.EnableFeasibilityMode = true;
% oc.fminopt.SubproblemAlgorithm = 'cg';
% oc.fminopt.TolProjCG = oc.fminopt.ConstraintTolerance;

% oc.optType = 'ipopt';
% oc.ipopt = struct;
% oc.ipopt.options = struct;
% oc.ipopt.options.ipopt = struct;
% oc.ipopt.options.ipopt.print_user_options = 'yes';
% oc.ipopt.options.ipopt.tol = 1e-4; % desired convergence tolerance
% oc.ipopt.options.ipopt.constr_viol_tol = 1e-4;
% oc.ipopt.options.ipopt.print_level = 5;
% oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
% oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';
% oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
% oc.ipopt.options.ipopt.mu_min = 1e-11;
% oc.ipopt.options.ipopt.mu_max = 1e-5;
% bound_push_frac = 1e-2;
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
[ ~, pulse, fields, opt ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, opt );

%% Prepare for optimization
opt.p0 = prev.opt.pOpt;
opt.pSc0 = prev.opt.pScOpt;

oc.pSc0 = opt.pSc0;
oc.p0 = opt.p0;

[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
if saveResult && strcmpi( solverName, "ipopt" )
    if ~isfolder( fullfile( oc.saveDir, testSubject ) )
        mkdir( fullfile( oc.saveDir, testSubject ) )
    end
    oc.ipopt.options.ipopt.output_file = char( fullfile( oc.saveDir, testSubject, "IPOPT_OUTPUT.txt" ) );
end
[ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

%% Run adjoint opt
clear costFnTrackDecisionVariablesWrapper;
clear costFnTrackConvergenceWrapper;
clear ipoptObjectiveBestOptWrapper;

opt = runAdjointOpt( opt, oc );

%% Post process optimization
oc.saveDir = fullfile( oc.saveDir, testSubject );
[ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
    opt, oc, pulse, fields );
