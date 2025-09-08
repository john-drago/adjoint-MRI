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
note = "Iterative universal nonselective pulse design - SPINS initialization";

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
    '2025-03-01-universal-nonselective',...
    '2025-03-01-universal-SPINS-init-90-cheb',...
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
st.useGPU = true;

st.binVis = false;
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
% numTrain = 13;
% numTest = 10;
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
pulse.universalTrain = sort( fields.universalTrain );
pulse.universalTest = sort( fields.universalTest );
% pulse.optPopulation = "tailored"; % e.g., universal or tailored

pulse.constantRotatingFrame = true;

% pulse.targFAVal = 10; % degrees
% pulse.flipAngleRangePlot = [ 0 20 ]; % Plot parameters
% pulse.dutyCycle = 0.05; % duty cycle in percent
% pulse.length = 250e-6; % time in ms
% st.opt_di = 0.005;
% st.opt_dt = 5e-6;
% st.val_di = 0.002;
% st.val_dt = 5e-6;

pulse.targFAVal = 90; % degrees
pulse.flipAngleRangePlot = [ 60 120 ]; % Plot parameters
pulse.dutyCycle = 0.025; % duty cycle in percent
pulse.length = 1e-3; % time in ms
st.opt_di = 0.005;
st.opt_dt = 10e-6;
st.val_di = 0.002;
st.val_dt = 10e-6;

pulse.targSt = struct;
pulse.targSt.x = fields.x;
pulse.targSt.y = fields.y;
pulse.targSt.z = fields.z;
pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
pulse.targSt.M0Vals = [ 0; 0; 1 ];

% Select below for the type of pulse

% SPINS specific parameters
pulse.min_kmax = 10; % rad/m
pulse.max_kmax = 100; % rad/m
pulse.min_u = ( 2*pi * 0.25 * 1000 ); % rad/m
pulse.max_u = ( 2*pi * 10 * 1000 ); % rad/m
pulse.min_v = ( 2*pi * 0.25 * 1000 ); % rad/m
pulse.max_v = ( 2*pi * 10 * 1000 ); % rad/m
pulse.min_a = 0.5;
pulse.max_a = 20;
pulse.min_b = 0;
pulse.max_b = 1;
pulse.initSlewTime = 30e-6;
pulse.finalSlewTime = 20e-6;

% % SPINS pulse
% pulse.name = "SPINS";
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

%% Initialize SPINS
initSPINSParams = struct;
initSPINSParams.kmax = 50; % rad / m
initSPINSParams.u = 8 * pi * 1000; % rad / s
initSPINSParams.v = 2 * pi * 1000; % rad / s
initSPINSParams.a = 10;
initSPINSParams.b = 0.5;

RFMax_constr_SPINSinit = 0.95 * optInit.RFMax_constr;
gradSlewRate_constr_SPINSinit = 0.95 * optInit.gradSlewRate_constr;

%% Initialize save struct
sp = struct;
sp.st = st;

sp.universalTrain = sort( universalTrain );
sp.universalTest = sort( universalTest );
sp.initOptTime_iters = zeros( numIter, 1 );

%% Solve SPINS
fprintf( "\n" );
fprintf( "-------------------------------------------------------------\n" );
fprintf( "Solving MLS SPINS initialization:\n" )
for ii = 1:numIter
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "iter:\t%i\n", ii );
    initTic = tic;
    [ RF0t_init, G0t_init, K0t_init, s0init, staStinit, outputinit ] = optimizeSPINSSTA(...
        initSPINSParams, RFMax_constr_SPINSinit, gradSlewRate_constr_SPINSinit,...
        optInit, pulse );
    optInit.initOptTime = toc( initTic );
    sp.initOptTime_iters( ii ) = optInit.initOptTime;
end
fprintf( "-------------------------------------------------------------\n" );
fprintf( "\n" );

RF0vec_init = reshape( transpose( RF0t_init ), [ optInit.numTimePoints * optInit.numXYCoils, 1 ] );

%% Generate SPINS waveform
optSPINS = struct;
optSPINS.structtype = 'opt';
optSPINS.di = oc.opt_di;
optSPINS.dt = round( oc.opt_dt, timeResDigits );
pulseSPINS = pulse;
pulseSPINS.name = "SPINS";
pulseSPINS.type = "base";
[ oc, pulseSPINS, fields, optSPINS ] = processFieldsStruct(...
    oc, pulseSPINS, fields, fields.si_train, optSPINS );

[ oc, pulseSPINS, optSPINS ] = processSPINSPulse_base( oc, pulseSPINS, optSPINS );
[ oc, pulseSPINS, optSPINS ] = processAdjointFunctions( oc, pulseSPINS, optSPINS );

optSPINS.p0 = zeros( optSPINS.numVars, 1 );
optSPINS.p0( optSPINS.breal_idx ) = real( RF0vec_init );
optSPINS.p0( optSPINS.bimag_idx ) = imag( RF0vec_init );
optSPINS.p0( optSPINS.kmax_idx ) = ( s0init.kmax - optSPINS.scVec( optSPINS.kmax_idx ) - optSPINS.min_kmax );
optSPINS.p0( optSPINS.a_idx ) = ( s0init.a - optSPINS.scVec( optSPINS.a_idx ) - optSPINS.min_a );
optSPINS.p0( optSPINS.b_idx ) = ( s0init.b - optSPINS.scVec( optSPINS.b_idx ) - optSPINS.min_b );
optSPINS.p0( optSPINS.u_idx ) = ( s0init.u - optSPINS.scVec( optSPINS.u_idx ) - optSPINS.min_u );
optSPINS.p0( optSPINS.v_idx ) = ( s0init.v - optSPINS.scVec( optSPINS.v_idx ) - optSPINS.min_v );

optSPINS.pSc0 = optSPINS.p0 ./ optSPINS.scVec;

wvSPINS = generateSPINSPlotWaveform_base( optSPINS.p0, optSPINS );

%% Project into SPINS basis
% optInit.p0 = zeros( optInit.numVars, 1 );
% optInit.p0( optInit.breal_idx ) = optSPINS.p0( optSPINS.breal_idx );
% optInit.p0( optInit.bimag_idx ) = optSPINS.p0( optSPINS.bimag_idx );
% optInit.p0( optInit.kmax_idx ) = ( s0init.kmax - optInit.scVec( optInit.kmax_idx ) - optInit.min_kmax );
% optInit.p0( optInit.a_idx ) = ( s0init.a - optInit.scVec( optInit.a_idx ) - optInit.min_a );
% optInit.p0( optInit.b_idx ) = ( s0init.b - optInit.scVec( optInit.b_idx ) - optInit.min_b );
% optInit.p0( optInit.u_idx ) = ( s0init.u - optInit.scVec( optInit.u_idx ) - optInit.min_u );
% optInit.p0( optInit.v_idx ) = ( s0init.v - optInit.scVec( optInit.v_idx ) - optInit.min_v );
% 
% optInit.pSc0 = optInit.p0 ./ optInit.scVec;
% oc.pSc0 = optInit.pSc0;
% oc.p0 = optInit.p0;

%% Project into PWC basis
% pwccoeffs = getPWCOptCoeffs( wvSPINS, optInit.tvec );
% 
% breal_idx_rshp = reshape( optInit.breal_idx, [ optInit.numTimePoints, optInit.numXYCoils ] );
% bimag_idx_rshp = reshape( optInit.bimag_idx, [ optInit.numTimePoints, optInit.numXYCoils ] );
% grad_idx_rshp = reshape( optInit.grad_idx, [ optInit.numTimePoints, 3 ] );
% 
% optInit.p0 = zeros( optInit.numVars, 1 );
% optInit.p0( breal_idx_rshp ) = pwccoeffs.breal;
% optInit.p0( bimag_idx_rshp ) = pwccoeffs.bimag;
% optInit.p0( grad_idx_rshp ) = pwccoeffs.grad;
% 
% optInit.pSc0 = optInit.p0 ./ optInit.scVec;
% 
% oc.pSc0 = optInit.pSc0;
% oc.p0 = optInit.p0;

%% Project into Cheb basis
chebTol = 1e-6;
chebcoeff = getChebOptCoeffs(...
    wvSPINS, chebTol );

pulse.numCheb_RF = min( size( chebcoeff.breal_coeffs, 1 ), optInit.orderCheb_RF + 1 );
pulse.numCheb_grad = min( size( chebcoeff.grad_coeffs, 1 ), optInit.orderCheb_grad + 1 );
pulse.numCheb_shim = 0;

pulse.orderCheb_RF = optInit.numCheb_RF - 1;
pulse.orderCheb_grad = optInit.numCheb_grad - 1;
pulse.orderCheb_shim = max( optInit.numCheb_shim - 1, 0 );

breal_idx_rshp = reshape( optInit.breal_idx, [ optInit.numCheb_RF, optInit.numXYCoils ] );
bimag_idx_rshp = reshape( optInit.bimag_idx, [ optInit.numCheb_RF, optInit.numXYCoils ] );
grad_idx_rshp = reshape( optInit.grad_idx, [ optInit.numCheb_grad, 3 ] );

p0 = zeros( optInit.numVars, 1 );
p0( breal_idx_rshp ) = chebcoeff.breal_coeffs( 1:pulse.numCheb_RF, : );
p0( bimag_idx_rshp ) = chebcoeff.bimag_coeffs( 1:pulse.numCheb_RF, : );
p0( grad_idx_rshp ) = chebcoeff.grad_coeffs( 1:pulse.numCheb_grad, : );

pSc0 = p0 ./ optInit.scVec;

optInit.p0 = p0;
optInit.pSc0 = pSc0;

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



            oc.fminopt

        elseif matches( solverName, 'ipopt' )

            oc.ipopt.options.ipopt.max_iter = maxOptIterationsIP;
            oc.ipopt.options.ipopt.max_wall_time = maxOptTime;



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

    clear oc opt;

end
