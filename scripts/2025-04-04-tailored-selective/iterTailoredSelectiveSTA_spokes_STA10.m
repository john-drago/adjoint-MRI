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
    "Iterative tailored slice-selective pulse design - spokes initialization";...
    "STA 10 degree flip angle, no further optimization" ];

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
    '2025-04-04-tailored-selective',...
    '2025-04-04-tailored-spokes-init-10-STA',...
    strcat( 'opt', timestr ) );

prevOptDirPath = fullfile( tld, "data", "opt",...
    "2025-04-03-universal-selective",...
    "2025-04-03-universal-spokes-init-10-pwc",...
    "XXXXXXXX",...
    "loc_-2.0cm",...
    "active-set" );

%% Initialize Script
st = struct;
st.numWorkers = 1;
st.binVis = false;
st.useGPU = true;
st.saveResult = true;
st.useParallel = false;
st.saveFileRecord = true;
st.ensureFeasibleStart = false;
st.trackConvergence = false;
st.trackDecisionVariables = false;
st.savePathTop = savePathTop;
st.currFile = currFile;
st.currDir = currDir;

st.numIter = 1;

st.val_di = 1e-3;
st.sliceThickness = 5e-3;
st.spokesFMinIterations = 300;
st.spokeLocations = 1e-2 * [...
    -2.0;...
    +0.0;...
    +2.0;...
    ];
st.numSpokeLocations = length( st.spokeLocations );

% 10 degrees
st.dt = 10e-6;
st.targFAVal = 10;
st.dutyCycle = 0.20;
st.flipAngleRangePlot = [ 0, 20 ];
st.numSpokes = 1;
st.centralSpokeTBW = 8;
st.centralSpokeLength = 1050 * 1e-6;
st.nonCentralSpokeTBW = 0;

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
    numTest = length( prev.valtest.subjIden );
    
    fprintf( "\n---------------------------------------------------------\n" );
    fprintf( "Time:\t%s\n", string( datetime ) );
    fprintf( "In directory:\t%s", folderPaths( dd ) );
    fprintf( "\n---------------------------------------------------------\n" );

    %% Iterate over slice locations
    for ll = 1:st.numSpokeLocations

        locString = sprintf( "loc_%+.1fcm", st.spokeLocations( ll ) * 1e2 );

        fprintf( "\n---------------------------------------------------------\n" );
        fprintf( "Time:\t%s\n", string( datetime ) );
        fprintf( "Location:\t%+.1fcm", st.spokeLocations( ll ) * 1e2 );
        fprintf( "\n---------------------------------------------------------\n" );

        %% Initialize save structs
        sp = struct;

        sp.optTime_iters = zeros( numTest, st.numIter );
        sp.fval_iters = zeros( numTest, st.numIter );

        sp.pOpt_iters = zeros( numTest, prev.opt.numVars );
        sp.pScOpt_iters = zeros( numTest, prev.opt.numVars );

        sp.subjIden = prev.valtest.subjIden;
        sp.optIden = optIden;
        sp.numIter = st.numIter;

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
            pulse.universalTrain = [];
            pulse.universalTest = [];
            
            pulse.dutyCycle = st.dutyCycle; % duty cycle in percent
            pulse.targFAVal = st.targFAVal; % degrees
            pulse.flipAngleRangePlot = st.flipAngleRangePlot;

            % spokes specific parameters
            pulse.sliceThickness = st.sliceThickness;
            pulse.sliceLocation = st.spokeLocations( ll );
            pulse.sliceDirection = [ 0; 0; 1 ];

            pulse.numSpokes = st.numSpokes;
            pulse.centralSpokeTBW = st.centralSpokeTBW;
            pulse.centralSpokeLength = st.centralSpokeLength;
            pulse.nonCentralSpokeTBW = st.nonCentralSpokeTBW;

            pulse.targSt = struct;
            pulse.targSt.M0Vals = [ 0; 0; 1 ];

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

            % PWC pulse
            pulse.name = "PWC";
            pulse.type = "base";

            %% Optimization Control (oc)
            st.saveDir = fullfile( savePathTop, locString, "STA" );

            [ oc, pulse ] = getOC2D( "active-set", pulse, fields, st );
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

            %% Run adjoint opt
            for ii = 1:st.numIter
                
                fprintf( "\n---------------------------------------------------------\n" );
                fprintf( "Time:\t%s\n", string( datetime ) );
                fprintf( "subject iter:\t%i", ii );
                fprintf( "\n---------------------------------------------------------\n" );

                %% Solve SPINS
                p = 8;
                roundTime = pulse.dt;
                initTic = tic;
                [ RF0init, K0init, ~, spokes ] = findInitialSpokes( opt, pulse,...
                    gradMaxConstr, gradSlewConstr, p, roundTime, st.spokesFMinIterations );

                [ RF0, G0, K0 ] = optimizeSpokesSTA(...
                    RF0init, K0init, RFMaxConstr, gradSlewConstr, opt, spokes, st.spokesFMinIterations );
                opt.initOptTime = toc( initTic );

                %% Generate spokes waveform
                opt.RF0 = RF0;
                opt.G0 = G0;
                opt.K0 = K0;

                wvSpokes = makeSpokesWaveformSample( RF0, G0, spokes );

                %% Project into PWC basis
                pwccoeffs = getPWCOptCoeffs( wvSpokes, opt.tvec );

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

                opt = runAdjointOpt( opt, oc );

                %% Track values
                sp.optTime_iters( tt, ii ) = opt.initOptTime;
                sp.fval_iters( tt, ii ) = opt.output.fval;

                sp.st = st;

                %% Post process optimization
                if ii == st.numIter
                    
                    sp.pOpt_iters( tt, : ) = opt.pOpt;
                    sp.pScOpt_iters( tt, : ) = opt.pScOpt;

                    oc.saveDir = fullfile( st.saveDir, sp.subjIden( tt ) );
                    
                    controlAdjointPostOpt( opt, oc, pulse, fields );
                    
                end
            end

            clear fields pulse opt oc valtrain valtest;

        end
        if st.saveResult
            save( fullfile( st.saveDir, "multipleSolutions.mat" ), '-struct', 'sp' );
        end

    end
end
