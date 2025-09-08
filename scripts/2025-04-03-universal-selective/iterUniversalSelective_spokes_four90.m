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
    "Iterative universal slice-selective pulse design - spokes initialization";...
    "Four-Cheb 90 degree flip angle" ];

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
    '2025-04-03-universal-selective',...
    '2025-04-03-universal-spokes-init-90-four',...
    strcat( 'opt', timestr ) );

st = struct;
st.currFile = currFile;
st.currDir = currDir;
st.savePathTop = savePathTop;
st.note = note;
st.numWorkers = 1;
st.useGPU = true;

st.binVis = false;
st.saveResult = true;
st.useParallel = false;
st.saveFileRecord = true;
st.ensureFeasibleStart = true;
st.trackConvergence = true;
st.trackDecisionVariables = true;
st.spacingTrackDecisionVariables = 25;

st.initialGuessOnly = false;

st.numIter = 1;
st.solverNames = [...
    "active-set";...
    "ipopt";...
    ];

st.val_di = 1e-3;
st.sliceThickness = 5e-3;
st.spokesFMinIterations = 300;

st.spokeLocations = 1e-2 * [...
    -2.0;...
    +0.0;...
    +2.0;...
    ];
st.numSpokeLocations = length( st.spokeLocations );

% 90 degrees
st.dt = 20e-6;
st.targFAVal = 90;
st.dutyCycle = 0.05;
st.flipAngleRangePlot = [ 0, 120 ];
st.numSpokes = 3;
st.centralSpokeTBW = 8;
st.centralSpokeLength = 1240 * 1e-6;
st.nonCentralSpokeTBW = 2;
st.RFBandwidth = 12.5e3;

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

%% Create optimization loop
for ll = 1:st.numSpokeLocations

    %% Fill-in fields struct
    fields = struct;
    fields.universalTrain = sort( universalTrain );
    fields.universalTest = sort( universalTest );
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

    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "Location:\t%+.1fcm", st.spokeLocations( ll ) * 1e2 );
    fprintf( "\n---------------------------------------------------------\n" );

    locString = sprintf( "loc_%+.1fcm", st.spokeLocations( ll ) * 1e2 );

    % ---------------------------------------------------------------------
    % Generate initial guess

    %% Define Pulse Structure
    pulse = struct;

    constraints = {...
        "total-RF-power", 100;... % W
        "max-RF-power", 24;... % W
        % "RF-slew-rate", 1e7;... % V/s
        % "RF-accel", 1.5e12;... % V/(sec^2)
        "RF-max", 2.190253975675163e+02;... % V
        % "RF-bandwidth", st.RFBandwidth;... % kHz
        % "peak-local-SAR", 1e4;... % W/kg
        % "peak-global-SAR", 1e4;... % W/kg
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
    pulse.optPopulation = "universal"; % e.g., universal or tailored
    pulse.universalTrain = sort( fields.universalTrain );
    pulse.universalTest = sort( fields.universalTest );
    % pulse.optPopulation = "tailored"; % e.g., universal or tailored

    pulse.dutyCycle = st.dutyCycle; % duty cycle in percent
    pulse.targFAVal = st.targFAVal; % degrees
    pulse.flipAngleRangePlot = st.flipAngleRangePlot;
    pulse.targSt = struct;
    pulse.targSt.M0Vals = [ 0; 0; 1 ];

    % spokes specific parameters
    pulse.sliceThickness = st.sliceThickness;
    pulse.sliceLocation = st.spokeLocations( ll );
    pulse.sliceDirection = [ 0; 0; 1 ];
    
    pulse.numSpokes = st.numSpokes;
    pulse.centralSpokeTBW = st.centralSpokeTBW;
    pulse.centralSpokeLength = st.centralSpokeLength;
    pulse.nonCentralSpokeTBW = st.nonCentralSpokeTBW;

    pulse.gyro = 267.5e6;

    RFMaxConstr = 0.95 * pulse.constraints( 'RF-max' );
    gradSlewConstr = 0.95 * pulse.constraints( 'grad-slew-rate' );
    gradMaxConstr = 0.95 * pulse.constraints( 'grad-max' );

    pulse.dt = st.dt;

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

    % Fourier-Cheb
    pulse.name = "fourier";
    pulse.type = "cheb";
    pulse.orderFourier_RF = 30;
    pulse.orderCheb_grad = 35;
    pulse.orderCheb_shim = 0;

    %% Optimization Control (oc)
    st.saveDir = savePathTop;
    [ oc, pulse ] = getOC2D( st.solverNames( 1 ), pulse, fields, st );

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

    %% Initialize save struct
    sp = struct;

    sp.universalTrain = sort( universalTrain );
    sp.universalTest = sort( universalTest ); 
    sp.initOptTime_iters = zeros( st.numIter, 1 );

    %% Initialize spokes
    p = 8;
    fprintf( "\n" );
    fprintf( "-------------------------------------------------------------\n" );
    fprintf( "Solving MLS Spokes initialization:\n" );
    for ii = 1:st.numIter
        fprintf( "iter:\t%i\n", ii );
        initTic = tic;
        [ RF0init, K0init, ~, spokes ] = findInitialSpokes( optInit, pulse,...
            gradMaxConstr, gradSlewConstr, p, roundTime, st.spokesFMinIterations );

        [ RF0, G0, K0 ] = optimizeSpokesSTA(...
            RF0init, K0init, RFMaxConstr, gradSlewConstr, optInit, spokes, st.spokesFMinIterations );
        optInit.initOptTime = toc( initTic );
        sp.initOptTime_iters( ii ) = optInit.initOptTime;
    end
    fprintf( "-------------------------------------------------------------\n" );
    fprintf( "\n" );

    %% Generate spokes waveform
    % RF0init_vec = reshape( RF0init, [ optInit.numXYCoils*pulse.numSpokes, 1 ] );
    % G0init_vec = reshape( G0init, [ 3*(pulse.numSpokes-1), 1 ] );

    % RF0_vec = reshape( RF0, [ optInit.numXYCoils*pulse.numSpokes, 1 ] );
    % G0_vec = reshape( G0, [ 3*(pulse.numSpokes-1), 1 ] );

    optInit.RF0 = RF0;
    optInit.G0 = G0;
    optInit.K0 = K0;
    optInit.spokes = spokes;

    wvSpokes = makeSpokesWaveformSample( RF0, G0, spokes );

    %% Project into Fourier-Cheb basis
    fouriercoeff = getFourierOptCoeffs(...
        wvSpokes, optInit.orderFourier_RF, [], [] );

    breal_idx_rshp = reshape( optInit.breal_idx, [ optInit.numFourier_RF, optInit.numXYCoils ] );
    bimag_idx_rshp = reshape( optInit.bimag_idx, [ optInit.numFourier_RF, optInit.numXYCoils ] );

    chebTol = 1e-6;
    chebcoeff = getChebOptCoeffs(...
        wvSpokes, chebTol );

    pulse.numCheb_grad = min( size( chebcoeff.grad_coeffs, 1 ), optInit.orderCheb_grad + 1 );
    pulse.numCheb_shim = 0;
    pulse.orderCheb_grad = optInit.numCheb_grad - 1;
    pulse.orderCheb_shim = max( optInit.numCheb_shim - 1, 0 );

    grad_idx_rshp = reshape( optInit.grad_idx, [ optInit.numCheb_grad, 3 ] );

    p0 = zeros( optInit.numVars, 1 );
    p0( breal_idx_rshp ) = fouriercoeff.breal_coeffs;
    p0( bimag_idx_rshp ) = fouriercoeff.bimag_coeffs;
    p0( grad_idx_rshp ) = chebcoeff.grad_coeffs( 1:pulse.numCheb_grad, : );

    pSc0 = p0 ./ optInit.scVec;

    optInit.p0 = p0;
    optInit.pSc0 = pSc0;

    oc.pSc0 = optInit.pSc0;
    oc.p0 = optInit.p0;

    %% Loop solvers

    for ss = 1:length( st.solverNames )

        %% Create opt struct
        opt = addOptInfoNewFieldStruct( optInit );

        %% Initialize OC
        if st.initialGuessOnly
            solverName = "STA";
            st.saveDir = fullfile( savePathTop, locString, solverName );
            
            st.trackConvergence = false;
            st.trackDecisionVariables = false;

            st.ensureFeasibleStart = false;

            [ oc, pulse ] = getOC2D( "active-set", pulse, fields, st );
            oc.fminopt.RelLineSrchBnd = 1e-13;
            oc.fminopt.MaxIterations = 0;
            oc.ensureFeasibleStart = false;

            % Initialize save structs
            sp.optTime_iters = zeros( 1, 1 );
            sp.fval_iters = zeros( 1, 1 );
            sp.numIter = st.numIter;

        else
            solverName = st.solverNames( ss );
            st.saveDir = fullfile( savePathTop, locString, solverName );
            [ oc, pulse ] = getOC2D( st.solverNames( ss ), pulse, fields, st );

            % oc.ensureFeasibleStartRelLineSrchBnd = 5e-5;
            
            maxOptTime = 30 * ( 60 * 60 );
            maxOptIterationsSQP = 2000;
            maxOptIterationsIP = 2000;

            if matches( solverName, 'active-set' )

                oc.fminopt.MaxIterations = maxOptIterationsSQP;
                oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );

                oc.fminopt.RelLineSrchBnd = 2.00e-4;
                
                oc.fminopt

            elseif matches( solverName, 'ipopt' )

                oc.ipopt.options.ipopt.max_iter = maxOptIterationsIP;
                oc.ipopt.options.ipopt.max_wall_time = maxOptTime;

                if abs( st.spokeLocations( ll ) - ( +2.0 * 1e-2 ) ) < 1e-12
                    
                    oc.constrTolSave = 3;

                    oc.ipopt.options.ipopt.max_iter = 2000;

                    oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
                    oc.ipopt.options.ipopt.mu_min = 1e-12;
                    oc.ipopt.options.ipopt.mu_max = 1e-7;
                    oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';

                    oc.ipopt.options.ipopt.tol = 1e-6;
                    oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
                    oc.ipopt.options.ipopt.dual_inf_tol = 1e-3;
                    oc.ipopt.options.ipopt.compl_inf_tol = 1e-3;

                    bound_push_frac = 1e-4;
                    oc.ipopt.options.ipopt.bound_push = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;

                    oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
                    oc.ipopt.options.ipopt.alpha_min_frac = 0.001;
                    oc.ipopt.options.ipopt.max_soc = 10;
                    oc.ipopt.options.ipopt.recalc_y = 'yes';
                    oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-6;

                    % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;

                    oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
                    oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
                    oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
                    oc.ipopt.options.ipopt.linear_solver = 'mumps';

                elseif abs( st.spokeLocations( ll ) - ( +0.0 * 1e-2 ) ) < 1e-12
                    
                    oc.constrTolSave = 3;

                    oc.ipopt.options.ipopt.max_iter = 2000;

                    oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
                    oc.ipopt.options.ipopt.mu_min = 1e-11;
                    oc.ipopt.options.ipopt.mu_max = 1e-5;
                    oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';

                    oc.ipopt.options.ipopt.tol = 1e-6;
                    oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
                    oc.ipopt.options.ipopt.dual_inf_tol = 1e-3;
                    oc.ipopt.options.ipopt.compl_inf_tol = 1e-3;

                    bound_push_frac = 1e-3;
                    oc.ipopt.options.ipopt.bound_push = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_relax_factor = 1e-3;

                    oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
                    oc.ipopt.options.ipopt.alpha_min_frac = 0.001;
                    oc.ipopt.options.ipopt.max_soc = 10;
                    oc.ipopt.options.ipopt.recalc_y = 'yes';
                    oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-6;

                    % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;

                    oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
                    oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
                    oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
                    oc.ipopt.options.ipopt.linear_solver = 'mumps';

                elseif abs( st.spokeLocations( ll ) - ( -2.0 * 1e-2 ) ) < 1e-12

                    oc.constrTolSave = 3;

                    oc.ipopt.options.ipopt.max_iter = 2000;

                    oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
                    oc.ipopt.options.ipopt.mu_min = 1e-12;
                    oc.ipopt.options.ipopt.mu_max = 1e-7;
                    oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';

                    oc.ipopt.options.ipopt.tol = 1e-6;
                    oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
                    oc.ipopt.options.ipopt.dual_inf_tol = 1e-3;
                    oc.ipopt.options.ipopt.compl_inf_tol = 1e-3;

                    bound_push_frac = 1e-4;
                    oc.ipopt.options.ipopt.bound_push = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
                    oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;

                    oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
                    oc.ipopt.options.ipopt.alpha_min_frac = 0.001;
                    oc.ipopt.options.ipopt.max_soc = 10;
                    oc.ipopt.options.ipopt.recalc_y = 'yes';
                    oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-6;

                    % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;

                    oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
                    oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
                    oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
                    oc.ipopt.options.ipopt.linear_solver = 'mumps'; 
                    
                end

                oc.ipopt.options.ipopt

            elseif matches( solverName, 'interior-point' )
                oc.fminopt.MaxIterations = maxOptIterationsIP;
                oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );

                oc.fminopt
            else
                error( "Unknown solver name." );
            end

            % Initialize save structs
            sp.numIter = st.numIter;

            sp.optTime_iters = zeros( st.numIter, 1 );
            sp.fval_iters = zeros( st.numIter, 1 );

            sp.track_fval_conv_iters = cell( st.numIter, 1 );
            sp.track_optTime_conv_iters = cell( st.numIter, 1 );
            sp.track_funccount_conv_iters = cell( st.numIter, 1 );
            sp.track_fvalraw_conv_iters = cell( st.numIter, 1 );
            sp.track_const_conv_iters = cell( st.numIter, 1 );

            sp.track_fval_vars_iters = cell( st.numIter, 1 );
            sp.track_optTime_vars_iters = cell( st.numIter, 1 );
            sp.track_funccount_vars_iters = cell( st.numIter, 1 );
            sp.track_pSc_vars_iters = cell( st.numIter, 1 );
            sp.track_p_vars_iters = cell( st.numIter, 1 );
            sp.track_fvalraw_vars_iters = cell( st.numIter, 1 );
            sp.track_const_vars_iters = cell( st.numIter, 1 );
        end

        oc.shimarray = true;

        oc.pSc0 = opt.pSc0;
        oc.p0 = opt.p0;

        fprintf( "\n---------------------------------------------------------\n" );
        fprintf( "Time:\t%s\n", string( datetime ) );
        fprintf( "Solver:\t%s", solverName );
        fprintf( "\n---------------------------------------------------------\n" );
        
        if ss > 1
            oc.ensureFeasibleStart = false;
        end
        
        [ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
        [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
        [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

        optInit.pSc0 = opt.pSc0;
        optInit.p0 = opt.p0;

        %% Save Structs

        for ii = 1:st.numIter

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

            %% Post process optimization

            if ( ii == st.numIter ) || st.initialGuessOnly
                
                controlAdjointPostOpt(... 
                    opt, oc, pulse, fields );

                if st.saveResult
                    sp.spokeLocation = st.spokeLocations( ll );
                    sp.st = st;
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

        clear oc opt;

    end

    clear pulse opt fields oc;

end