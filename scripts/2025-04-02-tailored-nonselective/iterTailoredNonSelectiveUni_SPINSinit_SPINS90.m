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
    "Iterative tailored nonselective pulse design - from universal pulse with SPINS initialization";...
    "SPINS 90 degree flip angle" ];

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

savePathTop = fullfile( tld, 'data', 'opt', ...
    '2025-04-02-tailored-nonselective',...
    '2025-04-02-tailored-SPINS-init-90-SPINS',...
    strcat( 'opt', timestr ) );

prevOptDirPath = fullfile( tld, "data", "opt", ...
    '2025-04-01-universal-nonselective',...
    "2025-04-01-universal-SPINS-init-90-SPINS",...
    "opt__250420043641" );

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

% st.opt_di = 0.015;
% st.opt_dt = 5e-6;
% st.val_di = 0.010;
% st.val_dt = 5e-6;

numIter = 1;
st.numIter = numIter;

%% Previous Optimization
[ ~, optIden ] = fileparts( prevOptDirPath );

% if multiple dirs to iterate over
prevOptDir = dir( prevOptDirPath );
prevOptDirNames = transpose( string( extractfield( prevOptDir, "name" ) ) );

if contains( prevOptDirNames, "optOutput.mat", 'IgnoreCase', true )
    folderPaths = prevOptDirPath;
elseif contains( prevOptDirNames, "180deg", 'IgnoreCase', true )
    folderPaths = fullfile( prevOptDirNames, "180deg" );
else
    prevOptDirLog = transpose( cell2mat( extractfield( prevOptDir, "isdir" ) ) );
    prevOptDir = prevOptDir( prevOptDirLog & ~matches( prevOptDirNames, [".", "..", ".DS_Store"], 'IgnoreCase', true ) );
    folderPaths = string(fullfile(prevOptDirPath, transpose( { prevOptDir.name } ) ) );
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

    if isfield( prev.oc, 'fminopt' )
        solverName = prev.oc.fminopt.Algorithm;
    elseif isfield( prev.oc, 'ipopt' )
        solverName = "ipopt";
    else
        error( "Unknown solver name." );
    end

    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "In directory:\t%s\n", folderPaths( dd ) );
    fprintf( "Using solver:\t%s", solverName );
    fprintf( "\n---------------------------------------------------------\n" );

    %% Initialize save structs
    sp = struct;
    sp.st = st;

    sp.optTime_iters = zeros( numTest, numIter );
    sp.fval_iters = zeros( numTest, numIter );

    sp.track_fval_conv_iters = cell( numTest, numIter );
    sp.track_optTime_conv_iters = cell( numTest, numIter );
    sp.track_funccount_conv_iters = cell( numTest, numIter );
    sp.track_fvalraw_conv_iters = cell( numTest, numIter );
    sp.track_const_conv_iters = cell( numTest, numIter );

    sp.track_fval_vars_iters = cell( numTest, numIter );
    sp.track_optTime_vars_iters = cell( numTest, numIter );
    sp.track_funccount_vars_iters = cell( numTest, numIter );
    sp.track_pSc_vars_iters = cell( numTest, numIter );
    sp.track_p_vars_iters = cell( numTest, numIter );
    sp.track_fvalraw_vars_iters = cell( numTest, numIter );
    sp.track_const_vars_iters = cell( numTest, numIter );

    sp.pOpt_iters = zeros( numTest, prev.opt.numVars );
    sp.pScOpt_iters = zeros( numTest, prev.opt.numVars );
    sp.p0_fromprev = prev.opt.pOpt;
    sp.pSc0_fromprev = prev.opt.pScOpt;
    
    sp.optIden = optIden;
    sp.subjIden = prev.valtest.subjIden;
    sp.universalTest = prev.valtest.subjIden;
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

        %% Optimization Control (oc)
        st.saveDir = ""; % will reassign shortly
        
        oc = getOC3D( solverName, st );

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
        st.saveDir = fullfile( savePathTop, solverName );

        oc = getOC3D( solverName, st );

        if tt > 1
            oc.ensureFeasibleStart = false;
        end

        maxOptTime = 1.0 * ( 60 );
        maxOptIterationsSQP = 1000;
        maxOptIterationsIP = 1000;

        if matches( solverName,  'active-set' )

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

        elseif matches( solverName, 'interior-point' )
            oc.fminopt.MaxIterations = maxOptIterationsIP;
            oc.fminopt.OutputFcn = prepareMaxTimeFMin( maxOptTime );

            oc.fminopt
        else
            error( "Unknown solver name." );
        end

        oc.shimarray = true;
        

        if tt == 1
            opt.p0 = prev.opt.pOpt;
            opt.pSc0 = prev.opt.pScOpt;
        else
            opt.p0 = p0;
            opt.pSc0 = pSc0;
        end

        oc.pSc0 = opt.pSc0;
        oc.p0 = opt.p0;

        [ oc, pulse, opt ] = processPulseName( oc, pulse, opt );
        [ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );

        if st.saveResult && strcmpi( solverName, "ipopt" )
            if ~isfolder( fullfile( st.saveDir, prev.valtest.subjIden( tt ) ) )
                mkdir( fullfile( st.saveDir, prev.valtest.subjIden( tt ) ) )
            end
            oc.ipopt.options.ipopt.output_file = char( fullfile( st.saveDir, prev.valtest.subjIden( tt ), "IPOPT_OUTPUT.txt" ) );
        end
        [ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

        if tt == 1
            p0 = opt.p0;
            pSc0 = opt.pSc0;
            sp.p0 = opt.p0;
            sp.pSc0 = opt.pSc0;
        end

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
            sp.optTime_iters( tt, ii ) = opt.optTime;
            sp.fval_iters( tt, ii ) = opt.output.fval;

            sp.track_fval_conv_iters{ tt, ii } = opt.fval_conv_iters;
            sp.track_optTime_conv_iters{ tt, ii } = opt.optTime_conv_iters;
            sp.track_funccount_conv_iters{ tt, ii } = opt.funccount_conv_iters;
            sp.track_fvalraw_conv_iters{ tt, ii } = opt.fvalraw_conv_iters;
            sp.track_const_conv_iters{ tt, ii } = opt.const_conv_iters;

            sp.track_fval_vars_iters{ tt, ii } = opt.fval_vars_iters;
            sp.track_optTime_vars_iters{ tt, ii } = opt.optTime_vars_iters;
            sp.track_funccount_vars_iters{ tt, ii } = opt.funccount_vars_iters;
            sp.track_pSc_vars_iters{ tt, ii } = opt.pSc_vars_iters;
            sp.track_p_vars_iters{ tt, ii } = opt.p_vars_iters;
            sp.track_fvalraw_vars_iters{ tt, ii } = opt.fvalraw_vars_iters;
            sp.track_const_vars_iters{ tt, ii } = opt.const_vars_iters;

            %% Post process optimization
            if ii == numIter
                sp.pOpt_iters( tt, : ) = opt.pOpt;
                sp.pScOpt_iters( tt, : ) = opt.pScOpt;

                oc.saveDir = fullfile( st.saveDir, prev.valtest.subjIden( tt ) );

                controlAdjointPostOpt( opt, oc, pulse, fields );

            end
        end

        if opt.useGPU
            gpuD = gpuDevice();
            reset( gpuD );
        end

        clear fields pulse opt oc valtrain valtest;

    end
    if st.saveResult
        save( fullfile( st.saveDir, "multipleSolutions.mat" ), '-struct', 'sp' );
    end
end
