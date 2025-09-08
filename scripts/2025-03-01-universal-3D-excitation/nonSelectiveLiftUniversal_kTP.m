%% File Initialization
% restoredefaultpath;
clear;
% close all;
home;

%% File Description
% Base script to optimize nonselective universal pulses using the adjoint
% method, while lifting to 180 degrees

%% Add Note To Describe Optimization
note = "Universal lifted nonselective pulse design - kTP initialization";

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
    '2025-03-01-universal-3D-excitation-lift-kTP-init-XXXName-XXXSolver',...
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

targFAVec = [...
    90;...
    120;...
    150;...
    180;...
    ];

numTargFAVec = length( targFAVec );

RFScaleVec = targFAVec / targFAVec( end );

finalRFMax = 2.190253975675163e+02;
finalAvgGlobalSAR = 3.2;
finalAvgLocalSAR = 20.0;
finalTotalRFPower = 100;
finalMaxRFPower = 24;

numIterOptInitial = 20;
numIterOptFinal = 30;

%% Name fields and path
fieldsPath = fullfile( tld, 'data', 'fields',...
    '2025-01-02-Siemens-Nova-8ch-pTx-Database-scaled', "UPdatabase.mat" );
VOPpath = fullfile( tld, 'data', 'fields',...
    '2025-01-02-Siemens-Nova-8ch-pTx-SAR', "Nova_8ch_pTx_SAR.mat" );

shimPath = [];

% load fields
UP = load( fieldsPath );
VOPst = load( VOPpath );

%% Get universal test and train data sets
% t = rng( "default" );
% numTrain = 5;
% numTest = 5;
% randSplit = randperm( length( UP.subjIden ) );
% 
% universalTrain = sort( UP.subjIden( randSplit( 1:numTrain ) ) );
% universalTest = sort( UP.subjIden( randSplit( (numTrain+1):(numTrain+numTest) ) ) );

universalTrain = sort( [...
    "AdjDataUser126";...
    "AdjDataUser109";...
    "AdjDataUser105";...
    "AdjDataUser120";...
    "AdjDataUser114";...
    "AdjDataUser110";...
    "AdjDataUser121";...
    "AdjDataUser117";...
    "AdjDataUser111";...
    "AdjDataUser108";...
    "AdjDataUser125";...
    "AdjDataUser123";...
    "AdjDataUser119";...
    ] );

universalTest = sort( [...
    "AdjDataUser103";...
    "AdjDataUser127";...
    "AdjDataUser104";...
    "AdjDataUser107";...
    "AdjDataUser122";...
    "AdjDataUser116";...
    "AdjDataUser112";...
    "AdjDataUser124";...
    "AdjDataUser113";...
    "AdjDataUser115";...
    ] );

idenLogical = contains( UP.subjIden, [ universalTrain; universalTest ], 'ignorecase', true );

% idenLogical = true;
fields = struct;
fields.universalTrain = sort( universalTrain );
fields.universalTest = sort( universalTest );

%% Fill-in fields struct
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
pulse = struct;

constraints = {...
    "total-RF-power", finalTotalRFPower * RFScaleVec( 1 )^2;... % W
    "max-RF-power", finalMaxRFPower * RFScaleVec( 1 )^2;... % W
    % "RF-slew-rate", 1e7;... % V/s
    % "RF-accel", 1.5e12;... % V/(sec^2)
    "RF-max", finalRFMax * RFScaleVec( 1 );... % V
    % "peak-local-SAR", 1e4;... % W/kg
    % "peak-global-SAR", 1e4;... % W/kg
    "average-local-SAR", finalAvgLocalSAR * RFScaleVec( 1 )^2;... % W/kg
    "average-global-SAR", finalAvgGlobalSAR * RFScaleVec( 1 )^2;... % W/kg
    "grad-max", 6e-3;... % T/m
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
pulse.excitationType = "non-selective";
pulse.terminalCostFunction = "arccos-least-squares";
pulse.runningCostFunction = [];

% specify tailored/universal identity if needed
pulse.optPopulation = "universal"; % e.g., universal or tailored
pulse.universalTrain = fields.universalTrain;
pulse.universalTest = fields.universalTest;
% pulse.optPopulation = "tailored"; % e.g., universal or tailored

pulse.constantRotatingFrame = true;

pulse.targFAVal = targFAVec( 1 ); % degrees
pulse.flipAngleRangePlot = [ 90 180 ]; % Plot parameters
pulse.targSt = struct;
pulse.targSt.x = fields.x;
pulse.targSt.y = fields.y;
pulse.targSt.z = fields.z;
pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
pulse.targSt.M0Vals = [ 0; 0; 1 ];

% Select below for the type of pulse

% kTP specific parameters
pulse.dutyCycle = 0.025; % duty cycle in percent
pulse.num_kTP = 9;
pulse.RFLength = 160e-6;
pulse.minRFSlewTime = 10e-6; % seconds, specific parameter for kTP pulses
pulse.blipLength = 70e-6;

pulse.length = pulse.num_kTP * pulse.RFLength + ( pulse.num_kTP - 1 ) * pulse.blipLength;

% % kTP pulse
% pulse.name = "kTP";
% pulse.type = "base";

% % PWC pulse
% pulse.name = "PWC"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
% pulse.type = "base"; % e.g., base or extended for kTP

% Cheb pulse
pulse.name = "cheb"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
pulse.type = "base"; % e.g., base or extended for kTP
pulse.orderCheb_RF = 20;
pulse.orderCheb_grad = 15;
pulse.orderCheb_shim = 0;

%% Optimization Control (oc)
oc = struct; % define opt control struct
oc.binVis = binVis; % whether or not to make images visible
oc.saveResult = saveResult; % whether or not to save optimization results
oc.saveDir = savePathTop; % location or where to save optimization results
oc.note = note;

oc.opt_di = 0.015 * ones( 3, 1 ); % dx, dy, dz for optimization
oc.opt_dt = 10e-6; % time in seconds

oc.val_di = 0.010 * ones( 3, 1 ); % dx, dy, dz for validation
oc.val_dt = 10e-6; % time in seconds

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
% oc.fminopt.TolConSQP = 1e-8;

% oc.fminopt.Algorithm = 'interior-point';
% oc.fminopt.StepTolerance = 1e-10;
% % oc.fminopt.BarrierParamUpdate = 'predictor-corrector';
% oc.fminopt.HessianApproximation = 'bfgs';
% oc.fminopt.InitBarrierParam = 1e-5;
% % oc.fminopt.EnableFeasibilityMode = true;
% oc.fminopt.SubproblemAlgorithm = 'cg';
% oc.fminopt.TolProjCG = oc.fminopt.ConstraintTolerance;

% oc.optType = "ipopt";
% oc.ipopt = struct;
% oc.ipopt.options = struct;
% oc.ipopt.options.ipopt = struct;
% oc.ipopt.options.ipopt.print_user_options = 'yes';
% % oc.ipopt.options.ipopt.print_options_documentation = 'yes';
% % oc.ipopt.options.ipopt.print_advanced_options = 'yes';
% oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
% oc.ipopt.options.ipopt.mu_min = 1e-11;
% oc.ipopt.options.ipopt.mu_max = 1e-5;
% oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';
% oc.ipopt.options.ipopt.tol = 1e-6;
% oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
% oc.ipopt.options.ipopt.dual_inf_tol = 1e-6;
% oc.ipopt.options.ipopt.compl_inf_tol = 1e-6;
% bound_push_frac = 1e-4;
% oc.ipopt.options.ipopt.bound_push = bound_push_frac;
% oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
% oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;
% oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
% oc.ipopt.options.ipopt.alpha_min_frac = 0.01;
% oc.ipopt.options.ipopt.max_soc = 10;
% oc.ipopt.options.ipopt.recalc_y = 'yes';
% oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-2;
% % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;
% oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
% oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
% oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
% oc.ipopt.options.ipopt.linear_solver = 'mumps';
% 
% oc.ipopt.options.ipopt.max_iter = maxOptIterations;
% oc.ipopt.options.ipopt.max_wall_time = maxOptTime

oc.saveFileRecord = saveFileRecord;
oc.currFile = currFile;
oc.useGPU = useGPU;
oc.ensureFeasibleStart = ensureFeasibleStart;
oc.trackConvergence = trackConvergence;
oc.trackDecisionVariables = trackDecisionVariables;
oc.spacingTrackDecisionVariables = 10;

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

%% Generate kTP opt struct
optkTP = struct;
optkTP.structtype = 'opt';
optkTP.di = oc.opt_di;
optkTP.dt = round( oc.opt_dt, timeResDigits );
pulsekTP = pulse;
pulsekTP.name = "kTP";
pulsekTP.type = "base";
[ oc, pulsekTP, fields, optkTP ] = processFieldsStruct(...
    oc, pulsekTP, fields, fields.si_train, optkTP );

[ oc, pulsekTP, optkTP ] = processkTPPulse_base( oc, pulsekTP, optkTP );
[ oc, pulsekTP, optkTP ] = processAdjointFunctions( oc, pulsekTP, optkTP );

%% Initialize kTP
gradSlewRate_constr_kTPinit = 0.95 * opt.gradSlewRate_constr;
RFMax_constr_kTPinit = 0.95 * opt.RFMax_constr;
p = 16;

initTic = tic;
% initial kTP
[ RF0init, G0init, K0init ] = findInitialkTPoint(...
    optkTP, pulse, RFMax_constr_kTPinit, gradSlewRate_constr_kTPinit, p );

RF0init_vec = reshape( RF0init, [ optkTP.numXYCoils*pulse.num_kTP, 1 ] );
G0init_vec = reshape( G0init, [ 3*(pulse.num_kTP-1), 1 ] );

% optimize kTP locations
[ RF0, G0, K0 ] = optimizekTPointsSTA(...
    RF0init, K0init, RFMax_constr_kTPinit, gradSlewRate_constr_kTPinit, optkTP, pulse );
opt.initOptTime = toc( initTic );

RF0_vec = reshape( RF0, [ optkTP.numXYCoils*pulse.num_kTP, 1 ] );
G0_vec = reshape( G0, [ 3*(pulse.num_kTP-1), 1 ] );

%% Generate kTP waveform
optkTP.p0 = zeros( optkTP.numVars, 1 );
optkTP.p0( optkTP.breal_idx ) = real( RF0_vec );
optkTP.p0( optkTP.bimag_idx ) = imag( RF0_vec );
optkTP.p0( optkTP.grad_idx ) = G0_vec;

optkTP.pSc0 = optkTP.p0 ./ optkTP.scVec;

wvkTP = generatekTPPlotWaveform_base( optkTP.p0, optkTP );

%% Project into kTP basis
% opt.p0 = zeros( opt.numVars, 1 );
% opt.p0( opt.breal_idx ) = optkTP.p0( optkTP.breal_idx );
% opt.p0( opt.bimag_idx ) = optkTP.p0( optkTP.bimag_idx );
% opt.p0( opt.grad_idx ) = optkTP.p0( optkTP.grad_idx );
% 
% opt.pSc0 = opt.p0 ./ opt.scVec;
% oc.pSc0 = opt.pSc0;
% oc.p0 = opt.p0;

%% Project into PWC basis
% pwccoeffs = getPWCOptCoeffs( wvkTP, opt.tvec );
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

%% Project into Cheb basis
chebTol = 1e-6;
chebcoeff = getChebOptCoeffs(...
    wvkTP, chebTol );

pulse.numCheb_RF = min( size( chebcoeff.breal_coeffs, 1 ), opt.orderCheb_RF + 1 );
pulse.numCheb_grad = min( size( chebcoeff.grad_coeffs, 1 ), opt.orderCheb_grad + 1 );
pulse.numCheb_shim = 0;

pulse.orderCheb_RF = opt.numCheb_RF - 1;
pulse.orderCheb_grad = opt.numCheb_grad - 1;
pulse.orderCheb_shim = max( opt.numCheb_shim - 1, 0 );

breal_idx_rshp = reshape( opt.breal_idx, [ opt.numCheb_RF, opt.numXYCoils ] );
bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numCheb_RF, opt.numXYCoils ] );
grad_idx_rshp = reshape( opt.grad_idx, [ opt.numCheb_grad, 3 ] );

p0 = zeros( opt.numVars, 1 );
p0( breal_idx_rshp ) = chebcoeff.breal_coeffs( 1:pulse.numCheb_RF, : );
p0( bimag_idx_rshp ) = chebcoeff.bimag_coeffs( 1:pulse.numCheb_RF, : );
p0( grad_idx_rshp ) = chebcoeff.grad_coeffs( 1:pulse.numCheb_grad, : );

pSc0 = p0 ./ opt.scVec;

opt.p0 = p0;
opt.pSc0 = pSc0;

oc.pSc0 = opt.pSc0;
oc.p0 = opt.p0;

%% Prepare for optimization
[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
[ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

%% Run adjoint opt

if trackConvergence
    opt.optTime_iters = zeros( numTargFAVec, 1 );
    opt.fval_iters = zeros( numTargFAVec, 1 );

    opt.track_fval_conv_iters = cell( numTargFAVec, 1 );
    opt.track_optTime_conv_iters = cell( numTargFAVec, 1 );
    opt.track_funccount_conv_iters = cell( numTargFAVec, 1 );
    
    opt.track_fval_vars_iters = cell( numTargFAVec, 1 );
    opt.track_optTime_vars_iters = cell( numTargFAVec, 1 );
    opt.track_funccount_vars_iters = cell( numTargFAVec, 1 );
    opt.track_pSc_vars_iters = cell( numTargFAVec, 1 );
    opt.track_p_vars_iters = cell( numTargFAVec, 1 );

    opt.pSc0_iters = zeros( numTargFAVec, opt.numVars );
    opt.p0_iters = zeros( numTargFAVec, opt.numVars );
    opt.pScOpt_iters = zeros( numTargFAVec, opt.numVars );
    opt.pOpt_iters = zeros( numTargFAVec, opt.numVars );
    opt.scVec_iters = zeros( numTargFAVec, opt.numVars );
    opt.targFAVec = targFAVec;
end

for ff = 1:numTargFAVec

    fprintf( "\ntarg FA:\t%g\n", targFAVec( ff ) );
    
    clear costFnTrackDecisionVariablesWrapper;
    clear costFnTrackConvergenceWrapper;
    clear ipoptObjectiveBestOptWrapper;

    if ff > 1
        pulse.targFAVal = targFAVec( ff );
        pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];

        Mtargvec = generateMtargOpt( pulse.targSt.MtargVals, fields.si_train, fields, opt );

        opt.Mtarg = Mtargvec;
        
        % scale RF
        opt.scVec( opt.breal_idx ) = finalRFMax / sqrt( 2 ) * RFScaleVec( ff );
        opt.scVec( opt.bimag_idx ) = finalRFMax / sqrt( 2 ) * RFScaleVec( ff );
        
        opt.p0( opt.breal_idx ) = opt.p0( opt.breal_idx ) * RFScaleVec( ff ) / RFScaleVec( ff - 1 );
        opt.p0( opt.bimag_idx ) = opt.p0( opt.bimag_idx ) * RFScaleVec( ff ) / RFScaleVec( ff - 1 );
        opt.pSc0 = opt.p0 ./ opt.scVec;

        oc.p0 = opt.p0;
        oc.pSc0 = opt.pSc0;

        % Scale Power
        opt.avgLocalSAR_constr = finalAvgLocalSAR * RFScaleVec( ff )^2;
        opt.avgGlobalSAR_constr = finalAvgGlobalSAR * RFScaleVec( ff )^2;
        opt.totalRFPower_constr = finalTotalRFPower * RFScaleVec( ff )^2;
        opt.maxRFPower_constr = finalMaxRFPower * RFScaleVec( ff )^2;
        opt.RFMax_constr = finalRFMax * RFScaleVec( ff );

        pulse.constraints( "average-local-SAR" ) = finalAvgLocalSAR * RFScaleVec( ff )^2;
        pulse.constraints( "average-global-SAR" ) = finalAvgGlobalSAR * RFScaleVec( ff )^2;
        pulse.constraints( "total-RF-power" ) = finalTotalRFPower * RFScaleVec( ff )^2;
        pulse.constraints( "max-RF-power" ) = finalMaxRFPower * RFScaleVec( ff )^2;
        pulse.constraints( "RF-max" ) = finalRFMax * RFScaleVec( ff );

    end

    if ff == numTargFAVec
        opt.fminopt.MaxIterations = numIterOptFinal;
    end
    
    % Prepare functions
    [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
    [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

    % run opt
    opt = runAdjointOpt( opt, oc );

    if trackConvergence
        opt.optTime_iters( ff ) = opt.optTime;
        opt.fval_iters( ff ) = opt.output.fval;
        opt.track_fval_conv_iters{ ff } = opt.fval_conv_iters;
        opt.track_optTime_conv_iters{ ff } = opt.optTime_conv_iters;
        opt.track_funccount_conv_iters{ ff } = opt.funccount_conv_iters;

        opt.track_fval_vars_iters{ ff } = opt.fval_vars_iters;
        opt.track_optTime_vars_iters{ ff } = opt.optTime_vars_iters;
        opt.track_funccount_vars_iters{ ff } = opt.funccount_vars_iters;
        opt.track_pSc_vars_iters{ ff } = opt.pSc_vars_iters;
        opt.track_p_vars_iters{  ff } = opt.p_vars_iters;

        opt.pSc0_iters( ff, : ) = opt.pSc0;
        opt.p0_iters( ff, : ) = opt.p0;
        opt.pScOpt_iters( ff, : ) = opt.pScOpt;
        opt.pOpt_iters( ff, : ) = opt.pOpt;
        opt.scVec_iters( ff, : ) = opt.scVec;
    end
    
    if ff < numTargFAVec
        opt.p0 = opt.pOpt;
        opt.pSc0 = opt.pScOpt;
    end

    if ff == numTargFAVec
        
        %% Post process optimization
        [ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
            opt, oc, pulse, fields );

    end

end


