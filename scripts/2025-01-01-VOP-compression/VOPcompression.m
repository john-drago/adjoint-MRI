%% File Initialization
restoredefaultpath;
clear;
close all;
home;

%% File Description
% This file will perform the VOP compression on the simulated fields from
% HFSS

%% Add paths to directories needed for optimization
currFile = strcat( mfilename( 'fullpath' ), ".m" );
currDir = fileparts( currFile );

cld = strsplit( currDir, 'scripts' );
cld = cld{ 1 };
addpath( genpath( fullfile( cld, 'util' ) ) );
addpath( genpath( fullfile( cld, 'adjoint' ) ) );
addpath( genpath( fullfile( cld, 'sim', 'SAR' ) ) );
tld = normalizePath( fullfile( cld, '..' ) );

dt = datetime('now');
dtIden = dt;
dtIden.Format = 'yyMMddHHmmss';
timestr = strcat('__', char( dtIden ) );

fprintf( "\n---------------------------------------------------------\n" );
fprintf( "Script start:\t%s\n", string(dt) );
fprintf( "Script Iden:\t%s", string(dtIden) );
fprintf( "\n---------------------------------------------------------\n" );

saveLoadDir = fullfile( tld, 'data', 'fields', ...
    '2025-01-01-simulated-fields' );

%% Change warning temporarily
toobigWarn = warning('error', 'MATLAB:save:sizeTooBigForMATFile');

%% Identify load files
loadFiles = [...
    "pTx_1x8_cylindrical_decoupled";...
    "BC2_back_feeding";...
    ];

%% Iterate over the load files
for ff = 1:length( loadFiles )
    
    %% Load File
    fileNameBase = loadFiles{ ff };
    em = load( fullfile( saveLoadDir, strcat( fileNameBase, ".mat" ) ) );
    
    %% Determine 10 g regions
    Ngrams = 10;
    [ NgramAssoc, NgramAmt, NgramDist ] = determineNgramRegions(...
        em.roi_body, em.rho_body, Ngrams, em.X, em.Y, em.Z );

    %% Compute Ngram Local SAR matrices
    [ QLocalMats, QLocalMatsIdx ] = determineNgramLocalSARMatrices(...
        NgramAssoc, em );

    %% Compute VOPs
    algorithm = 'kuehne';
    epsG = 0.03;
    overestMat = 'maxlocalSAReig';

    [ VOPs, VOPtoc, VOPID ] = performVOPCompression( QLocalMats, algorithm, epsG, overestMat );

    %% Compute Global SAR matrices
    [ QGlobalMat ] = determineGlobalSARMatrix( em );

    %% Save results
    ssVOP = struct;
    ssVOP.algorithm = algorithm;
    ssVOP.epsG = epsG;
    ssVOP.overestMat = overestMat;
    ssVOP.VOPs = VOPs;
    ssVOP.VOPtoc = VOPtoc;
    ssVOP.VOPID = VOPID;
    ssVOP.QGlobalMat = QGlobalMat;

    ssQLocal = struct;
    ssQLocal.QLocalMats = QLocalMats;
    ssQLocal.QLocalMatsIdx = QLocalMatsIdx;
    ssQLocal.QGlobalMat = QGlobalMat;
    ssQLocal.NgramAssoc = NgramAssoc;
    ssQLocal.NgramAmt = NgramAmt;
    ssQLocal.NgramDist = NgramDist;
    ssQLocal.Ngrams = Ngrams;

    %% Attempt save
    %% save VOPs
    saveDirVOP = fullfile( saveLoadDir, "VOPs" );
    if ~isfolder( saveDirVOP )
        mkdir( saveDirVOP );
    end

    try
        save( fullfile( saveDirVOP, strcat( fileNameBase, '.mat' ) ), '-struct',...
        'ssVOP' );
    catch
        save( fullfile( saveDirVOP, strcat( fileNameBase, '.mat' ) ), '-struct',...
        'ssVOP', "-v7.3" );
    end

    %% save QLocal
    saveDirQLocal = fullfile( saveLoadDir, "QLocal" );
    if ~isfolder( saveDirQLocal )
        mkdir( saveDirQLocal );
    end

    try
        save( fullfile( saveDirQLocal, strcat( fileNameBase, '.mat' ) ), '-struct',...
        'ssQLocal' );
    catch
        save( fullfile( saveDirQLocal, strcat( fileNameBase, '.mat' ) ), '-struct',...
        'ssQLocal', "-v7.3" );
    end

end

%% Change warning back
toobigWarn.state = 'on';
toobigWarn;
