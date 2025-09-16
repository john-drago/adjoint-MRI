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
note = [...
    "Iterative universal nonselective pulse design - kTP initialization";...
    "kTP 10 degree flip angle" ];

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

savePathTop = fullfile( tld, 'data', 'opt',...
    '2025-04-01-universal-nonselective',...
    '2025-04-01-universal-kTP-init-10-kTP',...
    strcat( 'opt', timestr ) );

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
st = struct;
st.numWorkers = 1;
st.binVis = false;
st.useGPU = true;
st.saveResult = true;
st.useParallel = false;
st.saveFileRecord = true;
st.ensureFeasibleStart = true;
st.trackConvergence = true;
st.trackDecisionVariables = true;
st.spacingTrackDecisionVariables = 25;
st.savePathTop = savePathTop;
st.currFile = currFile;
st.currDir = currDir;
st.note = note;

st.initialGuessOnly = false;

numIter = 1;
st.numIter = numIter;

solverNames = [...
    "active-set";...
    "ipopt";...
    ];

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
    "total-RF-power", 100;... % W
    "max-RF-power", 24;... % W
    % "RF-slew-rate", 1e7;... % V/s
    % "RF-accel", 1.5e12;... % V/(sec^2)
    "RF-max", 2.190253975675163e+02;... % V
    % "peak-local-SAR", 1e4;... % W/kg
    % "peak-global-SAR", 1e4;... % W/kg
    "average-local-SAR", 20;... % W/kg
    "average-global-SAR", 3.2;... % W/kg
    "grad-max", 4.95e-3;... % T/m
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
pulse.universalTrain = sort( fields.universalTrain );
pulse.universalTest = sort( fields.universalTest );
% pulse.optPopulation = "tailored"; % e.g., universal or tailored

pulse.constantRotatingFrame = true;

pulse.targFAVal = 10; % degrees
pulse.flipAngleRangePlot = [ 0 20 ]; % Plot parameters
pulse.dutyCycle = 0.05; % duty cycle in percent
pulse.num_kTP = 3;
pulse.RFLength = 50e-6;
pulse.minRFSlewTime = 10e-6; %% seconds, specific parameter for kTP pulses
pulse.blipLength = 50e-6;
st.opt_di = 0.005;
st.opt_dt = 5e-6;
st.val_di = 0.002;
st.val_dt = 5e-6;

pulse.targSt = struct;
pulse.targSt.x = fields.x;
pulse.targSt.y = fields.y;
pulse.targSt.z = fields.z;
pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
pulse.targSt.M0Vals = [ 0; 0; 1 ];

pulse.length = pulse.num_kTP * pulse.RFLength + ( pulse.num_kTP - 1 ) * pulse.blipLength;

% kTP pulse
pulse.name = "kTP";
pulse.type = "base";

%% Optimization Control (oc)
st.saveDir = savePathTop;
oc = getOC3D( solverNames( 1 ), st );

%% Initialize opt and val structs
optInit = struct;
optInit.structtype = 'opt';

timeResDigits = 6; % round to microseconds
optInit.di = oc.opt_di;
optInit.dt = round( oc.opt_dt, timeResDigits);

% Process fields struct 
% Get data about the fields such as number of coils and number of subjects
[ oc, pulse, fields ] = processFieldsData( oc, pulse, fields );

% Now add fields to the opt struct
[ oc, pulse, fields, optInit ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, optInit );

% Process Pulse struct
[ oc, pulse, optInit ] = processPulseName( oc, pulse, optInit );
[ oc, pulse, optInit ] = processAdjointFunctions( oc, pulse, optInit );

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

%% Initialize save struct
sp = struct;
sp.st = st;

sp.universalTrain = sort( universalTrain );
sp.universalTest = sort( universalTest );
sp.initOptTime_iters = zeros( numIter, 1 );

%% Initialize kTP
gradSlewRate_constr_kTPinit = 0.95 * optInit.gradSlewRate_constr;
RFMax_constr_kTPinit = 0.95 * optInit.RFMax_constr;
p = 16;

fprintf( "\n" );
fprintf( "-------------------------------------------------------------\n" );
fprintf( "Solving MLS kTP initialization:\n" )
for ii = 1:numIter
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "iter:\t%i\n", ii );
    initTic = tic;
    % initial kTP
    [ RF0init, G0init, K0init ] = findInitialkTPoint(...
        optkTP, pulse, RFMax_constr_kTPinit, gradSlewRate_constr_kTPinit, p );

    RF0init_vec = reshape( RF0init, [ optkTP.numXYCoils*pulse.num_kTP, 1 ] );
    G0init_vec = reshape( G0init, [ 3*(pulse.num_kTP-1), 1 ] );

    % optimize kTP locations
    [ RF0, G0, K0 ] = optimizekTPointsSTA(...
        RF0init, K0init, RFMax_constr_kTPinit, gradSlewRate_constr_kTPinit, optkTP, pulse );
    optInit.initOptTime = toc( initTic );
    sp.initOptTime_iters( ii ) = optInit.initOptTime;
end
fprintf( "-------------------------------------------------------------\n" );
fprintf( "\n" );

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
optInit.p0 = zeros( optInit.numVars, 1 );
optInit.p0( optInit.breal_idx ) = optkTP.p0( optkTP.breal_idx );
optInit.p0( optInit.bimag_idx ) = optkTP.p0( optkTP.bimag_idx );
optInit.p0( optInit.grad_idx ) = optkTP.p0( optkTP.grad_idx );

optInit.pSc0 = optInit.p0 ./ optInit.scVec;
oc.pSc0 = optInit.pSc0;
oc.p0 = optInit.p0;

%% Prepare for optimization
numSolvers = length( solverNames );

for ss = 1:length( solverNames )
    
    %% Create opt struct
    opt = addOptInfoNewFieldStruct( optInit );

    %% Initialize OC
    if st.initialGuessOnly

        solverName = "STA";
        st.saveDir = fullfile( savePathTop, solverName );

        st.trackConvergence = false;
        st.trackDecisionVariables = false;

        oc = getOC3D( "active-set", st );
        oc.fminopt.RelLineSrchBnd = 1e-13;
        oc.fminopt.MaxIterations = 0;

        oc.ensureFeasibleStart = false;
        
        % Initialize save structs
        sp.numIter = numIter;
        sp.optTime_iters = zeros( 1, 1 );
        sp.fval_iters = zeros( 1, 1 );

    else
        solverName = solverNames( ss );
        st.saveDir = fullfile( savePathTop, solverName );
        oc = getOC3D( solverName, st );

        if ss > 1
            oc.ensureFeasibleStart = false;
        end

        maxOptTime = 12 * ( 60*60 );
        maxOptIterationsSQP = 1500;
        maxOptIterationsIP = 750;

        if matches( solverName, "active-set" )

            oc.fminopt.MaxIterations = maxOptIterationsSQP;
            oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );

            oc.fminopt.RelLineSrchBnd = 5e-3;

            oc.fminopt

        elseif matches( solverName, 'ipopt' )

            oc.ipopt.options.ipopt.max_iter = maxOptIterationsIP;
            oc.ipopt.options.ipopt.max_wall_time = maxOptTime;

            oc.ipopt.options.ipopt.mu_min = 1e-11;
            oc.ipopt.options.ipopt.mu_max = 1e-5;

            oc.ipopt.options.ipopt.tol = 1e-6;
            oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
            oc.ipopt.options.ipopt.dual_inf_tol = 1e-6;
            oc.ipopt.options.ipopt.compl_inf_tol = 1e-6;

            bound_push_frac = 1e-4;
            oc.ipopt.options.ipopt.bound_push = bound_push_frac;
            oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
            oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;

            oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
            oc.ipopt.options.ipopt.alpha_min_frac = 0.01;
            oc.ipopt.options.ipopt.max_soc = 10;
            oc.ipopt.options.ipopt.recalc_y = 'yes';
            oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-2;

            % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;

            oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
            oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
            oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
            oc.ipopt.options.ipopt.linear_solver = 'mumps';

            oc.ipopt.options.ipopt

        elseif matches( solverName, "interior-point" )
            oc.fminopt.MaxIterations = maxOptIterationsIP;
            oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );

            oc.fminopt
        else
            error( "Unknown solver name." );
        end

        % Initialize save structs
        sp.numIter = numIter;
        sp.optTime_iters = zeros( numIter, 1 );
        sp.fval_iters = zeros( numIter, 1 );

        sp.track_fval_conv_iters = cell( numIter, 1 );
        sp.track_optTime_conv_iters = cell( numIter, 1 );
        sp.track_funccount_conv_iters = cell( numIter, 1 );
        sp.track_fvalraw_conv_iters = cell( numIter, 1 );
        sp.track_const_conv_iters = cell( numIter, 1 );

        sp.track_fval_vars_iters = cell( numIter, 1 );
        sp.track_optTime_vars_iters = cell( numIter, 1 );
        sp.track_funccount_vars_iters = cell( numIter, 1 );
        sp.track_pSc_vars_iters = cell( numIter, 1 );
        sp.track_p_vars_iters = cell( numIter, 1 );
        sp.track_fvalraw_vars_iters = cell( numIter, 1 );
        sp.track_const_vars_iters = cell( numIter, 1 );
    end

    oc.shimarray = true;
    

    oc.pSc0 = opt.pSc0;
    oc.p0 = opt.p0;

    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "Solver:\t%s", solverName );
    fprintf( "\n---------------------------------------------------------\n" );

    [ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
    [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
    [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

    optInit.pSc0 = opt.pSc0;
    optInit.p0 = opt.p0;

    %% Run adjoint opt
    for ii = 1:numIter

        fprintf( "\n---------------------------------------------------------\n" );
        fprintf( "Time:\t%s\n", string( datetime ) );
        fprintf( "iter:\t%i", ii );
        fprintf( "\n---------------------------------------------------------\n" );

        clear costFnTrackDecisionVariablesWrapper;
        clear costFnTrackConvergenceWrapper;
        clear ipoptObjectiveBestOptWrapper;

        opt = runAdjointOpt( opt, oc );

        %% Track values
        sp.optTime_iters( ii ) = opt.optTime; 
        sp.fval_iters( ii ) = opt.output.fval; 

        if ~st.initialGuessOnly
            sp.track_fval_conv_iters{ ii } = opt.fval_conv_iters;
            sp.track_optTime_conv_iters{ ii } = opt.optTime_conv_iters;
            sp.track_funccount_conv_iters{ ii } = opt.funccount_conv_iters;
            sp.track_fvalraw_conv_iters{ ii } = opt.fvalraw_conv_iters;
            sp.track_const_conv_iters{ ii } = opt.const_conv_iters;

            sp.track_fval_vars_iters{ ii } = opt.fval_vars_iters;
            sp.track_optTime_vars_iters{ ii } = opt.optTime_vars_iters;
            sp.track_funccount_vars_iters{ ii } = opt.funccount_vars_iters;
            sp.track_pSc_vars_iters{ ii } = opt.pSc_vars_iters;
            sp.track_p_vars_iters{ ii } = opt.p_vars_iters;
            sp.track_fvalraw_vars_iters{ ii } = opt.fvalraw_vars_iters;
            sp.track_const_vars_iters{ ii } = opt.const_vars_iters;
        else
            sp.optTime_iters( ii ) = 0;
        end

        sp.pSc0 = oc.pSc0;
        sp.p0 = oc.p0;
        sp.pScOpt = opt.pScOpt;
        sp.pOpt = opt.pOpt;

        sp.st = st;

        %% Post process optimization
        if ( ii == numIter ) || st.initialGuessOnly
            
            controlAdjointPostOpt( opt, oc, pulse, fields );

            if st.saveResult
                save( fullfile( st.saveDir, "multipleSolutions.mat" ), '-struct', 'sp' );
            end
        end

        if st.initialGuessOnly
                break;
        end
    end

    if st.initialGuessOnly
        break;
    end

    if opt.useGPU
        gpuD = gpuDevice();
        reset( gpuD );
        clear gpuD;
    end

end