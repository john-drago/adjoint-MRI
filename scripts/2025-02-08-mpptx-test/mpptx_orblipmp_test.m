%% File Initialization
restoredefaultpath;
clear;
close all;
home;

%% File Description
% This file will test the implementation procedure for the adjoint method
% (freely optimizable waveforms) with multiphoton parallel transmission
% paradigm (on-resonance, blip, and multiphoton subpulses).

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
timestr = strcat('__', char( dtIden ) );

fprintf( "\n---------------------------------------------------------\n" );
fprintf( "Script start:\t%s\n", string(dt) );
fprintf( "Script Iden:\t%s", string(dtIden) );
fprintf( "\n---------------------------------------------------------\n" );

savePathTop = fullfile( tld, 'data', 'opt', ...
    '2025-02-08-mpptx-orblipmp-test', strcat( 'opt', timestr ));

numWorkers = 1;
binVis = true;
useGPU = true;
saveResult = false;
useParallel = false;

%% Name fields and path
% subjIden = "pTx_1x8_cylindrical_decoupled";
subjIden = "BC2_back_feeding";
subjPath = fullfile( tld, 'data', 'fields', '2025-01-01-simulated-fields', subjIden );

fieldsPath = fullfile( strcat( subjPath, ".mat" ) );
VOPpath = fullfile( tld, 'data', 'fields', '2025-01-01-simulated-fields', 'VOPs',...
    strcat( subjIden, '.mat' ) );

shimPath = fullfile( tld, 'data', 'fields', '2025-01-03-sim-32ch-shim-coil',...
    '32ch_SensFields.mat');

% shimPath = [];

%% Initial Process of fields
UP = load( fieldsPath );
BS = load( shimPath );

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

%% Initialize fields struct
% fields struct should have the following attributes
%
%   - x discretization along x-axis (1D)
%   - y discretization along y-axis (1D)
%   - z discretization along z-axis (1D)
%   - X, Y, Z ndgrid format from x, y, z
%   - bzsens: fieldmaps from z-directed coils (not gradients) interpolated
%   to be at the points specified by X, Y, Z. The indices of bz will be:
%   [coil index, x index, y index, z index] if there are multiple coils and
%   [x index, y index, z index] if only one external z coil.
%   - b1p: xy-plane fieldmaps (complex-valued) interpolated at the points
%   specified by X, Y, Z. The indices will be: [coil index, x index, y 
%   index, z index, subj index] if there are multiple coils and [x index, y index, z 
%   index, subj index] if only one external z coil. If only one external
%   coil the indices could be [coil index, x index, y index, z index ].
%   - db0: the db0 value to use when performing optimization interpolated
%   to the points specified by X, Y, Z. If a shimmed db0 profile is
%   desired, then db0 should be the shimmed profile
%   - opt_roi: logical array of points upon which to perform optimization
%   interpolated at the points specified by X, Y, and Z
%   - val_roi: logical array of points upon which to perform validation
%   interpolated at the points specified by X, Y, and Z
%   - zi: this will be a Nsubj x 1 cell array that contains the z indices
%   (integers) for the roi points. this will help with plotting
%   - subjIden: this will be an Nsubj x 1 string array that contains
%   identifying information for the subject
%   - subjPath: this will be an Nsubj x 1 string array that contains the
%   path where each subject is stored

[ ~, ~, K ] = ndgrid( 1:length(UP.x), 1:length(UP.y), 1:length(UP.z) );

fields = struct;
fields.x = UP.x;
fields.y = UP.y;
fields.z = UP.z;
fields.X = UP.X;
fields.Y = UP.Y;
fields.Z = UP.Z;
fields.b1p = UP.b1p( :, :, :, idenLogical );
% fields.b1p = UP.b1p( :, :, :, :, idenLogical );
fields.db0 = UP.db0shim( :, :, :, idenLogical );
fields.opt_roi = UP.roi_brain( :, :, :, idenLogical );
fields.val_roi = UP.roi_body( :, :, :, idenLogical );

% fields.zi = UP.zi( idenLogical );
% fields.subjIden = string( UP.subjIden( idenLogical ) );
% fields.subjPath = string( UP.subjPath( idenLogical ) );

fields.zi = unique( K( UP.roi_brain ) );
fields.subjIden = subjIden;
fields.subjPath = subjPath;

fields.bz = zeros( [ BS.coilNum, size( fields.X ) ] );

for cc = 1:BS.coilNum
    fields.bz( cc, :, :, : ) = Interp3D(...
        BS.X, BS.Y, BS.Z, squeeze(BS.bz( cc, :, :, : )),...
        UP.X, UP.Y, UP.Z );
end
clear BS;
% fields.bz = 0;

% Save path names
fields.subjectMaps = fieldsPath;
fields.zCoilPath = shimPath;

%% Load VOPs
fields.VOPpath = VOPpath;
VOPst = load( VOPpath );
fields.VOPs = VOPst.VOPs;

% fields.QGlobal = VOPst.QPartialBody;
fields.QGlobal = VOPst.QGlobalMat;

%% Define pulse description
pulse = struct;

constraints = {...
    "peak-local-SAR", 1e4;... % W/kg
    "peak-global-SAR", 1e4;... % W/kg
    "average-local-SAR", 1e4;... % W/kg
    "average-global-SAR", 1e4;... % W/kg
    "total-RF-power", 1000;... % W
    "max-RF-power", 500;... % W
    "RF-accel", 5e12;... % V/(sec^2)
    "RF-slew-rate", 3e7;... % V/s
    "RF-max", 550;... % V
    "grad-max", 25e-3;... % T/m
    "grad-blip-max", 5e-3;...
    "grad-slew-rate", 200;... % T/m/sec
    "grad-accel", 1e7;... % T/m/sec^2
    "shim-max", 30;... % Amp-turns
    "shim-blip-max", 15;...
    "shim-total", 100;... % Amp-turns
    "shim-slew-rate", 3e6;... % Amp-turns/sec
    };

constraintName = string( constraints( :, 1 ) );
constraintValue = cell2mat( constraints( :, 2 ) );
pulse.constraints = dictionary( constraintName, constraintValue );

pulse.excitationType = "non-selective";
pulse.terminalCostFunction = "arccos-least-squares";
pulse.runningCostFunction = [];

% specify tailored/universal identity if needed
% pulse.optPopulation = "universal"; % e.g., universal or tailored
% pulse.universalTrain = universalTrain;
% pulse.universalTest = universalTest;

pulse.optPopulation = "tailored"; % e.g., universal or tailored
% pulse.tailoredTrain = tailoredTrain;

targFAVal = 10; % degrees
pulse.targSt = struct;
pulse.targSt.x = fields.x;
pulse.targSt.y = fields.y;
pulse.targSt.z = fields.z;
pulse.targSt.MtargVals = [ sind( targFAVal ); 0; cosd( targFAVal ) ];
pulse.targSt.M0Vals = [ 0; 0; 1 ];

pulse.dutyCycle = 0.05; % duty cycle

% mpptx specific parameters
pulse.name = "mpptx";
pulse.type = "or-blip-mp";

pulse.dfxy_mpptx = 5e3;
pulse.fz_mpptx = 5e3;
pulse.tORSP = 0.40e-3; % on-resonance subpulse time in seconds
pulse.tBlip = 0.10e-3; % on-resonance subpulse time in seconds
pulse.tMPSP = 0.50e-3; % on-resonance subpulse time in seconds
pulse.length = 1e-3; % time in ms
pulse.minSlewTime = 40e-6; % seconds, specific parameter for MP-pTx pulses

% Plot parameters
pulse.flipAngleRangePlot = [ 0 30 ];

oc = struct; % define opt control struct
oc.popt = struct;
oc.popt.zCoilScale = 0.4;
oc.popt.rfScale = 0.4;

%% Optimization Control (oc)
oc.binVis = binVis; % whether or not to make images visible
oc.saveResult = saveResult; % whether or not to save optimization results
oc.saveDir = savePathTop; % location or where to save optimization results

oc.opt_di = 0.015 * ones( 3, 1 ); % dx, dy, dz for optimization
oc.opt_dt = 20e-6; % time in seconds

oc.val_di = 0.010 * ones( 3, 1 ); % dx, dy, dz for validation
oc.val_dt = 20e-6; % time in seconds

% optimization specific parameters
oc.optType = "ga";

optDisplay = "iter";

oc.fminopt = optimoptions( "fmincon" );
oc.fminopt.Display = optDisplay;
oc.fminopt.UseParallel = useParallel;
oc.fminopt.MaxFunctionEvaluations = inf;
oc.fminopt.StepTolerance = 1e-4;
oc.fminopt.FiniteDifferenceStepSize = 1e-7;
oc.fminopt.FunctionTolerance = 1e-5; % doesn't matter for SQP
oc.fminopt.ConstraintTolerance = 1e-8;
oc.fminopt.SpecifyObjectiveGradient = true;
oc.fminopt.SpecifyConstraintGradient = true;
oc.fminopt.CheckGradients = false;

oc.fminopt.MaxIterations = 3;
oc.fminopt.Algorithm = 'active-set';
oc.fminopt.RelLineSrchBnd = 5e-3;
oc.fminopt.RelLineSrchBndDuration = inf;
oc.fminopt.TolConSQP = 1e-10;

oc.gaopt = optimoptions( "ga" );
oc.gaopt.Display = optDisplay;
oc.gaopt.UseParallel = useParallel;
oc.gaopt.MaxGenerations = 10;
oc.gaopt.MaxStallGenerations = 10;
oc.gaopt.PopulationSize = 200;
oc.gaopt.MaxTime = 60; % in seconds
oc.gaopt.MaxStallTime = 60; % in seconds
oc.gaopt.CrossoverFraction = 0.80;
oc.gaopt.MutationFcn = "mutationadaptfeasible";
oc.gaopt.CreationFcn = @gacreationnonlinearfeasible; % "gacreationlinearfeasible"
oc.gaopt.CrossoverFcn = "crossoverscattered";
oc.gaopt.SelectionFcn = "selectionroulette";
oc.gaopt.InitialPopulationRange = 0.005 * [ -1; 1 ];

if isfield( oc, 'rngstate')
    oc.gaopt.Type = oc.rngstate.Type;
    oc.gaopt.Seed = oc.rngstate.Seed;
    oc.gaopt.State = oc.rngstate.State;
end

oc.useGPU = useGPU;

%% Parallel Test
oc = determineParallelWorkers( numWorkers, oc );

%% Run Optimization
[ opt, valtrain, valtest, pulse, oc ] = controlAdjointOpt( oc, pulse, fields );


