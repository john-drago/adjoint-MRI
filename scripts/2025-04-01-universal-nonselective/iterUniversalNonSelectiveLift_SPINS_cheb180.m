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
    "Iterative universal nonselective 'lift' pulse design - SPINS initialization";...
    "Chebyshev 180 degree flip angle"];

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
    '2025-04-01-universal-SPINS-init-180-cheb',...
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

st.opt_di = 0.005;
st.opt_dt = 10e-6;
st.val_di = 0.002;
st.val_dt = 10e-6;

targFAVec = [...
    90;...
    100;...
    110;...
    120;...
    130;...
    140;...
    150;...
    160;...
    170;...
    175;...
    180;...
    ];

numTargFAVec = length( targFAVec );

RFScaleVec = targFAVec / targFAVec( end );

finalRFMax = 2.190253975675163e+02;
finalAvgGlobalSAR = 3.2;
finalAvgLocalSAR = 20.0;
finalTotalRFPower = 100;
finalMaxRFPower = 24;

numIterOptInitialIP = 500;
numIterOptMiddleIP = 250;
numIterOptFinalIP = 500;

numIterOptInitialSQP = 750;
numIterOptMiddleSQP = 300;
numIterOptFinalSQP = 750;
relLineSrchBndInitial = 1e-3;
relLineSrchBndFinal = 1e-3;

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
pulse.universalTrain = sort( fields.universalTrain );
pulse.universalTest = sort( fields.universalTest );
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

pulse.dutyCycle = 0.025; % duty cycle in percent
pulse.length = 2e-3; % time in ms

% Cheb pulse
pulse.name = "cheb"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
pulse.type = "base"; % e.g., base or extended for SPINS
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

%% Solve SPINS
fprintf( "\n" );
fprintf( "-------------------------------------------------------------\n" );
fprintf( "Solving MLS SPINS initialization:\n" )
fprintf( "-------------------------------------------------------------\n\n" );

initTic = tic;
[ RF0t_init, G0t_init, K0t_init, s0init, staStinit, outputinit ] = optimizeSPINSSTA(...
    initSPINSParams, RFMax_constr_SPINSinit, gradSlewRate_constr_SPINSinit,...
    optInit, pulse );
optInit.initOptTime = toc( initTic );

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
p0Init = optInit.p0;
pSc0Init = optInit.pSc0;

numSolvers = length( solverNames );

for ss = 1:length( solverNames )

    solverName = solverNames( ss );

    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "Solver:\t%s", solverName );
    fprintf( "\n---------------------------------------------------------\n" );

    st.saveDir = fullfile( savePathTop, solverName );
    
    %% Create opt struct
    opt = addOptInfoNewFieldStruct( optInit );

    %% Initialize OC
    oc = getOC3D( solverName, st );
    oc.shimarray = true;
    
    
    opt.pSc0 = pSc0Init;
    opt.p0 = p0Init;
    oc.pSc0 = pSc0Init;
    oc.p0 = p0Init;

    if ss > 1
        oc.ensureFeasibleStart = false;
    end

    [ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
    [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
    [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

    if ss == 1
        pSc0Init = opt.pSc0; 
        p0Init = opt.p0;
    end

    %% Initialize save structs
    sp = struct;
    sp.st = st;

    sp.universalTrain = sort( universalTrain );
    sp.universalTest = sort( universalTest );

    sp.initOptTime_iters = opt.initOptTime;

    sp.optTime_iters = zeros( numIter, numTargFAVec );
    sp.fval_iters = zeros( numIter, numTargFAVec );

    sp.track_fval_conv_iters = cell( numIter, numTargFAVec );
    sp.track_optTime_conv_iters = cell( numIter, numTargFAVec );
    sp.track_funccount_conv_iters = cell( numIter, numTargFAVec );
    sp.track_fvalraw_conv_iters = cell( numIter, numTargFAVec );
    sp.track_const_conv_iters = cell( numIter, numTargFAVec );

    sp.track_fval_vars_iters = cell( numIter, numTargFAVec );
    sp.track_optTime_vars_iters = cell( numIter, numTargFAVec );
    sp.track_funccount_vars_iters = cell( numIter, numTargFAVec );
    sp.track_pSc_vars_iters = cell( numIter, numTargFAVec );
    sp.track_p_vars_iters = cell( numIter, numTargFAVec );
    sp.track_fvalraw_vars_iters = cell( numIter, 1 );
    sp.track_const_vars_iters = cell( numIter, 1 );

    sp.pSc0_iters = zeros( numTargFAVec, opt.numVars );
    sp.p0_iters = zeros( numTargFAVec, opt.numVars );
    sp.pScOpt_iters = zeros( numTargFAVec, opt.numVars );
    sp.pOpt_iters = zeros( numTargFAVec, opt.numVars );
    sp.scVec_iters = zeros( numTargFAVec, opt.numVars );

    sp.targFAVec = targFAVec;
    sp.numTargFAVec = numTargFAVec;
    sp.numIter = numIter;

    sp.st = st;

    %% Run adjoint opt
    for ii = 1:numIter

        fprintf( "\n---------------------------------------------------------\n" );
        fprintf( "Time:\t%s\n", string( datetime ) );
        fprintf( "iter:\t%i", ii );
        fprintf( "\n---------------------------------------------------------\n" );

        for ff = 1:numTargFAVec

            fprintf( "\n---------------------------------------------------------\n" );
            fprintf( "Time:\t%s\n", string( datetime ) );
            fprintf( "Target FA:\t%g", targFAVec( ff ) );
            fprintf( "\n---------------------------------------------------------\n" );

            clear costFnTrackDecisionVariablesWrapper;
            clear costFnTrackConvergenceWrapper;
            clear ipoptObjectiveBestOptWrapper;

            oc.saveDir = fullfile( st.saveDir, sprintf( "%ideg",  targFAVec( ff ) ) );

            pulse.targFAVal = targFAVec( ff );
            pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
            Mtargvec = generateMtargOpt( pulse.targSt.MtargVals, fields.si_train, fields, opt );
            opt.Mtarg = Mtargvec;
            fields.Mtarg = zeros( [ size( fields.X ), 3 ] );
            fields.Mtarg( :, :, :, 1 ) = pulse.targSt.MtargVals( 1 );
            fields.Mtarg( :, :, :, 2 ) = pulse.targSt.MtargVals( 2 );
            fields.Mtarg( :, :, :, 3 ) = pulse.targSt.MtargVals( 3 );

            if ff == 1
                opt.pSc0 = pSc0Init;
                opt.p0 = p0Init;
                oc.pSc0 = pSc0Init;
                oc.p0 = p0Init;
            end

            % scale RF
            opt.scVec( opt.breal_idx ) = finalRFMax / sqrt( 2 ) * RFScaleVec( ff );
            opt.scVec( opt.bimag_idx ) = finalRFMax / sqrt( 2 ) * RFScaleVec( ff );

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

            if ( ff > 1 )

                opt.p0( opt.breal_idx ) = opt.p0( opt.breal_idx ) * RFScaleVec( ff ) / RFScaleVec( ff - 1 );
                opt.p0( opt.bimag_idx ) = opt.p0( opt.bimag_idx ) * RFScaleVec( ff ) / RFScaleVec( ff - 1 );
                opt.pSc0 = opt.p0 ./ opt.scVec;

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

            oc.p0 = opt.p0;
            oc.pSc0 = opt.pSc0;

            maxOptTime = 12 * ( 60*60 );

            if matches( solverName, 'active-set' )

                if ff == numTargFAVec
                    oc.fminopt.MaxIterations = numIterOptFinalSQP;
                elseif ff == 1
                    oc.fminopt.MaxIterations = numIterOptInitialSQP;
                else
                    oc.fminopt.MaxIterations = numIterOptMiddleSQP;
                end
                if ( ff == numTargFAVec ) || ( ff == ( numTargFAVec - 1 ) )
                    oc.fminopt.RelLineSrchBnd = relLineSrchBndFinal;
                else
                    oc.fminopt.RelLineSrchBnd = relLineSrchBndInitial;
                end

                oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );



                oc.fminopt

            elseif matches( solverName, 'ipopt' )

                if ff == numTargFAVec
                    oc.ipopt.options.ipopt.max_iter = numIterOptFinalIP;
                elseif ff == 1
                    oc.ipopt.options.ipopt.max_iter = numIterOptInitialIP;
                else
                    oc.ipopt.options.ipopt.max_iter = numIterOptMiddleIP;
                end
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

            elseif matches( solverName, 'interior-point' )

                if ff == numTargFAVec
                    oc.fminopt.MaxIterations = numIterOptFinalIP;
                elseif ff == 1
                    oc.fminopt.MaxIterations = numIterOptInitialIP;
                else
                    oc.fminopt.MaxIterations = numIterOptMiddleIP;
                end
                oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );


                oc.fminopt
            else
                error( "Unknown solver name." );
            end

            if ( ff == 1 )
                oc.ensureFeasibleStart = false;
            else
                oc.ensureFeasibleStart = true;
            end

            % Prepare functions
            [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );
            [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

            opt = runAdjointOpt( opt, oc );

            %% Track values
            sp.optTime_iters( ii, ff ) = opt.optTime;
            sp.fval_iters( ii, ff ) = opt.output.fval;

            sp.track_fval_conv_iters{ ii, ff } = opt.fval_conv_iters;
            sp.track_optTime_conv_iters{ ii, ff } = opt.optTime_conv_iters;
            sp.track_funccount_conv_iters{ ii, ff } = opt.funccount_conv_iters;
            sp.track_fvalraw_conv_iters{ ii, ff } = opt.fvalraw_conv_iters;
            sp.track_const_conv_iters{ ii, ff } = opt.const_conv_iters;

            sp.track_fval_vars_iters{ ii, ff } = opt.fval_vars_iters;
            sp.track_optTime_vars_iters{ ii, ff } = opt.optTime_vars_iters;
            sp.track_funccount_vars_iters{ ii, ff } = opt.funccount_vars_iters;
            sp.track_pSc_vars_iters{ ii, ff } = opt.pSc_vars_iters;
            sp.track_p_vars_iters{ ii, ff } = opt.p_vars_iters;
            sp.track_fvalraw_vars_iters{ ii, ff } = opt.fvalraw_vars_iters;
            sp.track_const_vars_iters{ ii, ff } = opt.const_vars_iters;

            if ( ii == numIter )
                sp.pSc0_iters( ff, : ) = opt.pSc0;
                sp.p0_iters( ff, : ) = opt.p0;
                sp.pScOpt_iters( ff, : ) = opt.pScOpt;
                sp.pOpt_iters( ff, : ) = opt.pOpt;
                sp.scVec_iters( ff, : ) = opt.scVec;
            end

            if ( ii == numIter )

                %% Post process optimization
                controlAdjointPostOpt( opt, oc, pulse, fields );

            end

            if ff < numTargFAVec
                opt.p0 = opt.pOpt;
                opt.pSc0 = opt.pScOpt;
            end

            if ( ff == numTargFAVec ) && ( ii == numIter )

                if st.saveResult
                    save( fullfile( st.saveDir, "multipleSolutions.mat" ), '-struct', 'sp' );
                end

            end

        end

    end

    if opt.useGPU
        gpuD = gpuDevice();
        reset( gpuD );
        clear gpuD;
    end


end