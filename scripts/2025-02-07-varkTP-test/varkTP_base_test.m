%% File Initialization
restoredefaultpath;
clear;
close all;
home;

%% File Description
% This file will test the implementation procedure for the adjoint method
% (freely optimizable waveforms) with variable kTP basis (blips).

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
    '2025-02-07-varkTP-base-test', strcat( 'opt', timestr ));

numWorkers = 1;
binVis = true;
useGPU = true;
saveResult = false;
useParallel = false;

%% Name fields and path
subjIden = "pTx_1x8_cylindrical_decoupled";
% subjIden = "BC2_back_feeding";
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
% fields.b1p = UP.b1p( :, :, :, idenLogical );
fields.b1p = UP.b1p( :, :, :, :, idenLogical );
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
    "max-pulse-length", 2e-3;... % seconds
    "total-RF-power", 1e3;... % W
    "max-RF-power", 1e3;... % W
    "RF-slew-rate", 5.5e7;... % V/s
    % "RF-accel", 5e12;... % V/(sec^2)
    "RF-max", 550;... % V
    "peak-local-SAR", 10e3;... % W/kg
    "peak-global-SAR", 10e3;... % W/kg
    "average-local-SAR", 10;... % W/kg
    "average-global-SAR", 3.2;... % W/kg
    "grad-max", 6e-3;... % T/m
    "grad-slew-rate", 200;... % T/m/sec
    "shim-max", 30;... % Amp-turns
    "shim-total", 150;... % Amp-turns
    "shim-slew-rate", 3e6;... % Amp-turns/sec
    };

constraintName = string( constraints( :, 1 ) );
constraintValue = cell2mat( constraints( :, 2 ) ); 
pulse.constraints = dictionary( constraintName, constraintValue );

pulse.optPopulation = "tailored"; % e.g., universal or tailored
pulse.excitationType = "non-selective"; % e.g., non-selective, slice-selective, inversion, refocusing,...
% pulse.terminalCostFunction = "magnitude-least-squares"; % e.g., magntiude-least-squares, mz-flip-angle, mxy-flip-angle, least-squares, weighted-magnitude-least-squares,...
pulse.terminalCostFunction = "arccos-least-squares";
pulse.runningCostFunction = [];

targFAVal = 10; % degrees
pulse.targSt = struct;
pulse.targSt.x = fields.x;
pulse.targSt.y = fields.y;
pulse.targSt.z = fields.z;
pulse.targSt.MtargVals = [ sind( targFAVal ); 0; cosd( targFAVal ) ];
pulse.targSt.M0Vals = [ 0; 0; 1 ];

% Plot parameters
pulse.flipAngleRangePlot = [ 0 20 ];

pulse.dutyCycle = 0.05; % duty cycle in percent

% var ktp specific parameters
pulse.name = "variable-ktp";
pulse.type = "base";

pulse.num_kTP = 5;

pulse.initRFLength = 120e-6;
pulse.initBlipLength = 60e-6;

pulse.minRFSlewTime = 10e-6; % seconds, specific parameter for kTP pulses

pulse.minGradSlewTime = 5e-6; % seconds, specific parameter for kTP pulses
pulse.minShimSlewTime = 5e-6; % seconds, specific parameter for kTP pulses

pulse.length = pulse.num_kTP * pulse.initRFLength + ...
    ( pulse.num_kTP - 1 ) * pulse.initBlipLength;


%% Optimization Control (oc)
oc = struct; % define opt control struct
oc.binVis = binVis; % whether or not to make images visible
oc.saveResult = saveResult; % whether or not to save optimization results
oc.saveDir = savePathTop; % location or where to save optimization results

oc.opt_di = 0.015 * ones( 3, 1 ); % dx, dy, dz for optimization
oc.opt_dt = 5e-6; % time in seconds

oc.val_di = 0.010 * ones( 3, 1 ); % dx, dy, dz for validation
oc.val_dt = 5e-6; % time in seconds

% optimization specific parameters
oc.optType = "fmin";

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
oc.gaopt.MaxGenerations = 2;
oc.gaopt.MaxStallGenerations = 2;
oc.gaopt.PopulationSize = 50;
oc.gaopt.MaxTime = inf; % in seconds
oc.gaopt.MaxStallTime = inf; % in seconds
oc.gaopt.CrossoverFraction = 0.80;
oc.gaopt.MutationFcn = "mutationadaptfeasible";
oc.gaopt.CreationFcn = @gacreationnonlinearfeasible; % "gacreationlinearfeasible"
oc.gaopt.CrossoverFcn = "crossoverscattered";
oc.gaopt.SelectionFcn = "selectionroulette";
oc.gaopt.InitialPopulationRange = 0.025 * [ -1; 1 ];

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


