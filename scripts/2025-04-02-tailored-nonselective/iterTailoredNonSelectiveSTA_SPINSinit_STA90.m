%% File Initialization
% restoredefaultpath;
clear;
% close all;
home;

%% File Description
% Base script to optimize nonselective tailored pulses using the adjoint
% method
%
% Specific pulse type:

%% Add Note To Describe Optimization
note = [...
    "Iterative tailored nonselective pulse design - SPINS initialization";...
    "STA 90 degree flip angle, no further optimization" ];

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
    '2025-04-02-tailored-nonselective',...
    '2025-04-02-tailored-SPINS-init-90-STA',...
    strcat( 'opt', timestr ) );

prevOptDirPath = fullfile( tld, "data", "opt",...
    "2025-04-01-universal-nonselective",...
    "2025-04-01-universal-SPINS-init-90-pwc",...
    "XXXXXXXX",...
    "active-set" );

%% Initialize Script
st = struct;
st.numWorkers = 1;
st.binVis = false;
st.useGPU = true;
st.saveResult = true;
st.useParallel = false;
st.saveFileRecord = true;
st.trackConvergence = false;
st.trackDecisionVariables = false;
st.savePathTop = savePathTop;
st.ensureFeasibleStart = false;
st.currFile = currFile;
st.currDir = currDir;

% st.opt_di = 0.015;
% st.opt_dt = 5e-6;
% st.val_di = 0.010;
% st.val_dt = 5e-6;

numIter = 1;

%% Previous Optimization
[ ~, optIden ] = fileparts( prevOptDirPath );

% if multiple dirs to iterate over
prevOptDir = dir( prevOptDirPath );
prevOptDirNames = transpose( string( extractfield( prevOptDir, "name" ) ) );

if ~contains( prevOptDirNames, "optOutput.mat", 'IgnoreCase', true )
    prevOptDirLog = transpose( cell2mat( extractfield( prevOptDir, "isdir" ) ) );
    prevOptDir = prevOptDir( prevOptDirLog & ~matches( prevOptDirNames, [ ".", "..", ".DS_Store" ], 'IgnoreCase', true ) );
    folderPaths = string(fullfile(prevOptDirPath, transpose( { prevOptDir.name } ) ) );
else
    folderPaths = prevOptDirPath;
end

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
for dd = 1:length( folderPaths )

    prev = load( fullfile( folderPaths( dd ), 'optOutput.mat' ) );

    st.opt_di = prev.oc.opt_di;
    st.opt_dt = prev.oc.opt_dt;
    st.val_di = prev.oc.val_di;
    st.val_dt = prev.oc.val_dt;

    numTest = length( prev.valtest.subjIden );

    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "In directory:\t%s", folderPaths( dd ) );
    fprintf( "\n---------------------------------------------------------\n" );

    %% Initialize save structs
    sp = struct;
    sp.st = st;

    sp.optTime_iters = zeros( numTest, numIter );
    sp.fval_iters = zeros( numTest, numIter );

    sp.pOpt_iters = zeros( numTest, prev.opt.numVars );
    sp.pScOpt_iters = zeros( numTest, prev.opt.numVars );

    sp.subjIden = prev.valtest.subjIden;
    sp.optIden = optIden;
    sp.numIter = numIter;

    %% Iterate over the different valtest subjects
    for tt = 1:numTest

        fprintf( "\n---------------------------------------------------------\n" );
        fprintf( "Time:\t%s\n", string( datetime ) );
        fprintf( "Test subject:\t%s", prev.valtest.subjIden( tt ) );
        fprintf( "\n---------------------------------------------------------\n" );

        %% Get validation data set
        idenLogical = contains( UP.subjIden, prev.valtest.subjIden( tt ), 'ignorecase', true );

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

        % PWC pulse
        pulse.name = "PWC"; % e.g., mpptx, instant, pseudo-spectral-shim, pseudo-spectral,...
        pulse.type = "base"; % e.g., base or extended for kTP

        %% Optimization Control (oc)
        st.saveDir = fullfile( savePathTop, "STA" );

        oc = getOC3D( "active-set", st );
        oc.fminopt.RelLineSrchBnd = 1e-13;
        oc.fminopt.MaxIterations = 0;
        oc.ensureFeasibleStart = false;

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

        % Process pulse struct
        [ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
        [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );

        %% Initialize SPINS
        initSPINSParams = struct;
        initSPINSParams.kmax = 50; % rad / m
        initSPINSParams.u = 8 * pi * 1000; % rad / s
        initSPINSParams.v = 2 * pi * 1000; % rad / s
        initSPINSParams.a = 10;
        initSPINSParams.b = 0.5;

        RFMax_constr_SPINSinit = 0.95 * opt.RFMax_constr;
        gradSlewRate_constr_SPINSinit = 0.95 * opt.gradSlewRate_constr;

        %% Run adjoint opt
        for ii = 1:numIter

            fprintf( "\n---------------------------------------------------------\n" );
            fprintf( "Time:\t%s\n", string( datetime ) );
            fprintf( "iter:\t%i", ii );
            fprintf( "\n---------------------------------------------------------\n" );

            %% Solve SPINS
            initTic = tic;
            [ RF0t_init, G0t_init, K0t_init, s0init, staStinit, outputinit ] = optimizeSPINSSTA(...
                initSPINSParams, RFMax_constr_SPINSinit, gradSlewRate_constr_SPINSinit,...
                opt, pulse );
            opt.initOptTime = toc( initTic );

            RF0vec_init = reshape( transpose( RF0t_init ), [ opt.numTimePoints * opt.numXYCoils, 1 ] );

            %% Generate SPINS waveform
            optSPINS = struct;
            optSPINS.structtype = 'opt';
            optSPINS.di = oc.opt_di;
            optSPINS.dt = round( oc.opt_dt, timeResDigits );
            pulseSPINS = pulse;
            pulseSPINS.name = "SPINS";
            pulseSPINS.type = "base";
            [ ~, pulseSPINS, fields, optSPINS ] = processFieldsStruct(...
                oc, pulseSPINS, fields, fields.si_train, optSPINS );

            [ ~, pulseSPINS, optSPINS ] = processSPINSPulse_base( oc, pulseSPINS, optSPINS );
            [ ~, pulseSPINS, optSPINS ] = processAdjointFunctions( oc, pulseSPINS, optSPINS );

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

            %% Project into PWC basis
            pwccoeffs = getPWCOptCoeffs( wvSPINS, opt.tvec );

            breal_idx_rshp = reshape( opt.breal_idx, [ opt.numTimePoints, opt.numXYCoils ] );
            bimag_idx_rshp = reshape( opt.bimag_idx, [ opt.numTimePoints, opt.numXYCoils ] );
            grad_idx_rshp = reshape( opt.grad_idx, [ opt.numTimePoints, 3 ] );

            opt.p0 = zeros( opt.numVars, 1 );
            opt.p0( breal_idx_rshp ) = pwccoeffs.breal;
            opt.p0( bimag_idx_rshp ) = pwccoeffs.bimag;
            opt.p0( grad_idx_rshp ) = pwccoeffs.grad;

            opt.pSc0 = opt.p0 ./ opt.scVec;

            oc.pSc0 = opt.pSc0;
            oc.p0 = opt.p0;

            %% Prepare for optimization
            [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

            clear costFnTrackDecisionVariablesWrapper;
            clear costFnTrackConvergenceWrapper;
            clear ipoptObjectiveBestOptWrapper;

            opt = runAdjointOpt( opt, oc );

            %% Track values
            sp.optTime_iters( tt, ii ) = opt.initOptTime;
            sp.fval_iters( tt, ii ) = opt.output.fval;

            %% Post process optimization
            if ii == numIter
                sp.pOpt_iters( tt, : ) = opt.pOpt;
                sp.pScOpt_iters( tt, : ) = opt.pScOpt;

                oc.saveDir = fullfile( st.saveDir, prev.valtest.subjIden( tt ) );
                [ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
                    opt, oc, pulse, fields );
            end
        end

    end
    if st.saveResult
        save( fullfile( st.saveDir, "multipleSolutions.mat" ), '-struct', 'sp' );
    end
end
