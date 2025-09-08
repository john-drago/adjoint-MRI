%% File Initialization
restoredefaultpath;
clearvars;
% close all;
home;

%% Description
% This script will load in the adjust data fieldmaps provided by Siemens
% and process them into a database that can used for pulse optimization
% 
% Script must be run on a linux computer than has BET, freesurfer, gzip
% installed
%% Files included
%
% Script was prepared to process the following adjust users: 
% AdjDataUser103.mat
% AdjDataUser104.mat
% AdjDataUser105.mat
% AdjDataUser107.mat
% AdjDataUser108.mat
% AdjDataUser109.mat
% AdjDataUser110.mat
% AdjDataUser111.mat
% AdjDataUser112.mat
% AdjDataUser113.mat
% AdjDataUser114.mat
% AdjDataUser115.mat
% AdjDataUser116.mat
% AdjDataUser117.mat
% xxxxxxxxxxxxxxxxxx AdjDataUser118.mat
% AdjDataUser119.mat
% AdjDataUser120.mat
% AdjDataUser121.mat
% AdjDataUser122.mat
% AdjDataUser123.mat
% AdjDataUser124.mat
% AdjDataUser125.mat
% AdjDataUser126.mat
% AdjDataUser127.mat

%% Set Paths 
currFile = strcat( mfilename( 'fullpath' ), ".m" );
currDir = fileparts( currFile );

cld = strsplit( currDir, 'scripts' );
cld = cld{ 1 };
addpath( genpath( fullfile( currDir, 'func' ) ) );
addpath( genpath( fullfile( cld, 'util' ) ) );
addpath( genpath( fullfile( cld, 'sim', 'b0shim' ) ) );
tld = normalizePath( fullfile( cld, '..' ) );

dt = datetime('now');
dtIden = dt;
dtIden.Format = 'yyMMddHHmmss';
timestr = strcat('__', char( dtIden ) );

fprintf( "\n---------------------------------------------------------\n" );
fprintf( "Script start:\t%s\n", string(dt) );
fprintf( "Script Iden:\t%s", string(dtIden) );
fprintf( "\n---------------------------------------------------------\n" );

%% File control
binVis = false;
saveFigs = true;
scaleb1p = true;
gyro = 267.5e6; % rad / ( s * T )
brainDens = 1040; % kg / m^3

%% Set save dir
saveDirPath = fullfile(...
    tld,...
    'data', 'fields',...
    '2025-01-02-Siemens-Nova-8ch-pTx-Database-scaled' );
saveName = "UPdatabase.mat";

%% Change warning temporarily
toobigWarn = warning('error', 'MATLAB:save:sizeTooBigForMATFile');

%% List Directory for loading
loadDir = fullfile( saveDirPath, "AdjExamples" );
lDir = dir( fullfile( loadDir, "*.mat" ) );
lDirNames = string( { lDir( : ).name }.' );
lDir = lDir( ~contains( lDirNames, "AdjDataUser118", "IgnoreCase", true ) );

numFM = length( lDir );

%% Define interpolation domain
% This is the domain in which we want to represent the pulses
% DCS coordinate system:
% x: right -> left ( - -> + )
% y: posterior -> anterior ( - -> + )
% z: superior -> inferior ( - -> + )
%
% Assuming that B0 points in the +z direction

dxn = 0.002;
xn = -0.160 : dxn : 0.160;
yn = -0.160 : dxn : 0.160;
zn = -0.200 : dxn : 0.200;

n_xn = length( xn );
n_yn = length( yn );
n_zn = length( zn );

[ Xn, Yn, Zn ] = ndgrid( xn, yn, zn );
Rn = cat( 4, Xn, Yn, Zn );

%% Initialize mean tracker
meanb1pSubj = zeros( numFM, 1 );

%% Iterate across the files
for ff = 1:numFM

    %% Load file
    ls = load( fullfile( lDir( ff ).folder, strcat( lDir( ff ).name ) ) );

    if isfield( ls, 'Adj' )
        Adj = ls.Adj;
    else
        Adj = ls;
    end

    if ff == 1
        b1p = zeros( Adj.coils, n_xn, n_yn, n_zn, numFM );
        b1p_CP = zeros( n_xn, n_yn, n_zn, numFM );
        db0 = zeros( n_xn, n_yn, n_zn, numFM );
        db0shim = zeros( n_xn, n_yn, n_zn, numFM );
        roi_body = false( n_xn, n_yn, n_zn, numFM );
        roi_brain = false( n_xn, n_yn, n_zn, numFM );
        zi = cell( numFM, 1 );
        weightMask = zeros( n_xn, n_yn, n_zn, numFM );
        subjIden = strings( numFM, 1 );
        subjPath = strings( numFM, 1 );
    end

    %% Get fields
    b1p_mns = reshape( Adj.S, [ Adj.image_m, Adj.image_n, Adj.coils, Adj.slices ] );
    b1p_mns = permute( b1p_mns, [ 3, 1, 2, 4 ] );

    db0_mns = Adj.B0;

    bodyMask_old = Adj.W;

    %% Get coordinates in (phase-encode)-(read-out)-(slice-select) frame
    % pe_pos = Adj.values_m / 1000; % convert to meters
    % ro_pos = Adj.values_n / 1000; % convert to meters
    % ss_pos = Adj.values_s / 1000; % convert to meters
    % [ PE, RO, SS ] = ndgrid( pe_pos, ro_pos, ss_pos );

    % %% Get positions ndgrid
    [ M, N, S ] = ndgrid( Adj.values_m, Adj.values_n, Adj.values_s );
    M = 1e-3 * M; % convert from mm to m; M is PE
    N = 1e-3 * N; % convert from mm to m; N is RO
    S = 1e-3 * S; % convert from mm to m; S is SS

    %% Coordinate System Handling
    % We are trying to get all coordinates in the DCS coordinate
    % system:
    %
    % DCS coordinate system:
    % x: right -> left ( - -> + )
    % y: posterior -> anterior ( - -> + )
    % z: superior -> inferior ( - -> + )
    %
    % Assuming that B0 points in the +z direction
    %
    % In the following code, we define a ( rotation / reflection ) matrix
    % that will be used to project the values in the coordinate system of
    % the image to the DCS coordinate system.
    %
    % In other words, the ( rotation / reflection ) matrix describes the
    % positions of the *NEW* axis (DCS) in the coordinate system of the
    % *OLD* axis (axis of the image).

    if Adj.image_ori == 0 % axial orientation
        % assume that the axial orientation is as follows:
        % M: phase-encode: Anterior (-) --> Posterior (+): y-axis
        % N: read-out: Left (-) --> Right (+): x-axis
        % S: slice-select: Inferior (-) --> Superior (+): z-axis
        %
        %
        R = [ [ 0; -1; 0], [ -1; 0; 0], [ -1; 0; 0] ];

        % Assume that the x-y axes in the axial orientation (as described
        % above) must become the x-y axes in the DCS, respectively. The x-y axes in the
        % axial orientation are rotated 180 degrees CCW about the -z-axis
        % in the axial orientation to become the x-y axes in the DCS.
        b1ptrans_fn = @( b1p ) b1p .* exp( 1j * -( 180 * pi/180 ) );

    elseif Adj.image_ori == 1 % sagittal orientation
        % assume that the sagittal orientation is as follows:
        % M: phase-encode: Anterior (-) --> Posterior (+): x-axis
        % N: read-out: Inferior (-) --> Superior (+): y-axis
        % S: slice-select: Right (-) --> Left (+): z-axis
        %
        %
        R = [ [ 0; 0; 1], [ -1; 0; 0], [ 0; -1; 0] ];
    
        % Assume that the z-y axes in the sagittal orientation (as
        % described above) must become the x-y axes in the DCS,
        % respectively. The y-z axes in the sagittal orientation are
        % reflected about the z axis in the sagittal orientation to become
        % the x-y axes in the DCS. y-axis in sagittal orientation stays the
        % y-axis in the DCS and the z-axis in the coronal orientation
        % becomes the x-axis in the DCS.
        b1ptrans_fn = @( b1p ) conj( b1p ) .* exp( 1j *  2 * ( 0 * pi/180 ) );

    elseif Adj.image_ori == 2 % coronal orientation
        % assume that the coronal orientation is as follows:
        % M: phase-encode: Inferior (-) --> Superior (+): x-axis
        % N: read-out: Right (-) --> Left (+): y-axis
        % S: slice-select: Anterior (-) --> Posterior (+): z-axis
        %
        %
        R = [ [ 0; 1; 0], [ 0; 0; -1], [ -1; 0; 0] ];

        % Assume that the x-z axes in the coronal orientation (as described
        % above) must become the x-y axes in the DCS, respectively. The x-z
        % axes in the coronal orientation are reflected about the x-axis in
        % the coronal orientation to become the x-y axes in the DCS. x-axis
        % in coronal orientation becomes the x-axis in the DCS and the
        % z-axis in the coronal orientation becomes the y-axis in the DCS
        % (after reflection).
        b1ptrans_fn = @( b1p ) conj( b1p ) .* exp( 1j *  2 * ( 0 * pi/180 ) );

    else
        error( "Unknown image orientation." )
    end

    %% Project locations in new coordinate system to position in the old
    % coordinate system
    r_newpos_old = transpose( R * transpose( [ Xn(:), Yn(:), Zn(:) ] ) );

    %% interpolate values at new locations in old coordinate system non-vector
    % brainMask_newpos_old = interpScalar( M, N, S, double( brainMask_old ), r_newpos_old );
    bodyMask_newpos_old = interpScalar( M, N, S, double( bodyMask_old ), r_newpos_old );
    db0_newpos_old = interpScalar( M, N, S, double( bodyMask_old ) .* double( db0_mns ), r_newpos_old );
    % weight_newpos_old = interpScalar( M, N, S, double( weightMask ), r_newpos_old );

    %% Get indices in new coordinate system
    in = 1 : length(xn);
    jn = 1 : length(yn);
    kn = 1 : length(zn);

    [ In, Jn, Kn ] = ndgrid( in, jn, kn );

    ip = In(:);
    jp = Jn(:);
    kp = Kn(:);

    scind = sub2ind( size(Xn), ip, jp, kp );

    %% Put scalar vectors into ndgrid array in the new coordinate system
    bodyMask_new = false( size( Xn ) );
    % brainMask_new = false( size( Xn ) );
    db0_iter = zeros( size( Xn ) );
    % weightMask_iter = zeros( size( Xn ) );

    logThresh = 0.5;
    bodyMask_new( scind ) = bodyMask_newpos_old > logThresh;
    % brainMask_new( scind ) = brainMask_newpos_old > logThresh;
    % weightMask_iter( scind ) = weight_newpos_old;
    
    % if any( ( R( :, 3 ) ) < 0 ) % corresponds to bottom row in inverse matrix
    %     sc = -1;
    % else
    %     sc = 1;
    % end
    sc = 1;
    
    %% Make db0 map
    db0_iter( scind ) = ( ( sc * 2 * pi ) / gyro ) * db0_newpos_old; % scale from Hz to T

    %% Smooth the db0 map
    db0_iter = smoothKernel( db0_iter );

    %% interpolate projection values at new locations in old coordinate system vector
    numCoils = size( b1p_mns, 1 );
    b1p_iter = zeros( [ numCoils, size( Xn ) ] );

    for nn = 1:numCoils

        b1p_newpos_old = interpScalar( M, N, S, double( bodyMask_old ) .* squeeze( b1p_mns(nn,:,:,:) ), r_newpos_old );

        %% Project vector values into the new coordinate system
        b1p_proj_unsrt = b1ptrans_fn( b1p_newpos_old );

        %% Get indices in new coordinate system
        ptxscind = sub2ind( [ numCoils, size(Xn) ], nn*ones(size(ip)), ip, jp, kp );

        %% Put b1p into ndgrid array in the new coordinate system
        b1p_iter( ptxscind ) = b1p_proj_unsrt * 1e-6; % scale from uT/V to T/V
        b1p_real_iter = smoothKernel( squeeze( real( b1p_iter( nn, :, :, : ) ) ) );
        b1p_imag_iter = smoothKernel( squeeze( imag( b1p_iter( nn, :, :, : ) ) ) );
        b1p_iter( nn, :, :, : ) = complex( b1p_real_iter, b1p_imag_iter );

    end

    %% Make CP mode with Nova coil
    CPangles = linspace( 0, 2*pi, ( Adj.coils + 1) ).';
    CPanglesComp = exp( 1j * CPangles( 1:(end-1) ) );
    CPanglesFOV_b1p = permute(...
        repmat( reshape( CPanglesComp, [ 1, 1, 1, Adj.coils ] ), size( b1p_iter, [ 2, 3, 4 ] ) ),...
        [ 4, 1, 2, 3 ] );
    b1p_CP_iter = squeeze( sum( b1p_iter .* CPanglesFOV_b1p, 1 ) );

    %% Generate brain mask with b0 using ROMEO idea
    [ brainMask_new_unsmooth, weightMask_iter ] = generateMaskFromB0(...
        db0_iter, double( bodyMask_new ), b1p_CP_iter );

    %% Erode and Dilate brain mask
    numPix = 1;
    dilPix = 9;
    [ brainMask_new ] = erodeBrainROI( brainMask_new_unsmooth, numPix, dilPix );

    %% Ensure body mask contains brain mask
    bodyMask_new = bodyMask_new | brainMask_new;
    
    %% Find the center of mass in the coordinate system
    mass_iter = double( brainMask_new ) .* ( brainDens * dxn^3 );
    r_com = squeeze( sum( Rn .* mass_iter, [ 1, 2, 3 ] ) )/ ( sum( mass_iter, "all" ) );

    dzi_com = round( r_com( 3 ) / dxn );

    %% Shift pixels according to center of mass to make isocenter at the center of mass
    zi_body = unique( Kn( bodyMask_new ) );
    zi_body = min( zi_body ) : max( zi_body );
    zi_body_shift = zi_body - dzi_com;

    % initialize all arrays
    bodyMask_new_shift = false( size( bodyMask_new ) );
    brainMask_new_shift = false( size( brainMask_new ) );
    b1p_iter_shift = zeros( size( b1p_iter ) );
    b1p_CP_iter_shift = zeros( size( b1p_CP_iter ) );
    db0_iter_shift = zeros( size( db0_iter ) );
    weightMask_iter_shift = zeros( size( weightMask_iter ) );

    % add the shifted data
    bodyMask_new_shift( :, :, zi_body_shift ) = bodyMask_new( :, :, zi_body );
    brainMask_new_shift( :, :, zi_body_shift ) = brainMask_new( :, :, zi_body );
    b1p_iter_shift( :, :, :, zi_body_shift ) = b1p_iter( :, :, :, zi_body );
    b1p_CP_iter_shift( :, :, zi_body_shift ) = b1p_CP_iter( :, :, zi_body );
    db0_iter_shift( :, :, zi_body_shift ) = db0_iter( :, :, zi_body );
    weightMask_iter_shift( :, :, zi_body_shift ) = weightMask_iter( :, :, zi_body );
    
    % get rid of shift arrays
    bodyMask_new = bodyMask_new_shift;
    brainMask_new = brainMask_new_shift;
    b1p_iter = b1p_iter_shift;
    b1p_CP_iter = b1p_CP_iter_shift;
    db0_iter = db0_iter_shift;
    weightMask_iter = weightMask_iter_shift;

    clear bodyMask_new_shift brainMask_new_shift b1p_iter_shift b1p_CP_iter_shift db0_iter_shift weightMask_iter_shift;

    %% Shim db0
    db0shim_iter = calcSH2DB0Correction( db0_iter, brainMask_new, bodyMask_new, Xn, Yn, Zn );

    %% Determine the z indices of interest
    zi_opt_buffer = 5;
    zi_opt_iter = unique( Kn( brainMask_new ) );
    zi_val_iter = unique( Kn( bodyMask_new ) );
    zi_opt_exp = max( [ ( min( zi_opt_iter ) - zi_opt_buffer ), min( zi_val_iter ) ] ) : min( [ ( max( zi_opt_iter ) + zi_opt_buffer ), max( zi_val_iter ) ] );

    %% Adjust the masks
    bodyMask_iter = false( size( bodyMask_new ) );
    brainMask_iter = false( size( brainMask_new ) );

    bodyMask_iter( :, :, zi_opt_exp ) = bodyMask_new( :, :, zi_opt_exp );
    brainMask_iter( :, :, zi_opt_exp ) = brainMask_new( :, :, zi_opt_exp );

    clear bodyMask_new brainMask_new;

    %% Determine mean abs( b1p_CP )
    meanb1pSubj( ff ) = mean( abs( b1p_CP_iter( brainMask_iter ) ) );

    %% Assign to arrays
    b1p( :, :, :, :, ff ) = b1p_iter;
    b1p_CP( :, :, :, ff ) = b1p_CP_iter;
    db0( :, :, :, ff ) = db0_iter;
    db0shim( :, :, :, ff ) = db0shim_iter;
    roi_body( :, :, :, ff ) = bodyMask_iter;
    roi_brain( :, :, :, ff ) = brainMask_iter;
    weightMask( :, :, :, ff ) = weightMask_iter;
    zi{ ff } = zi_opt_iter;
    subjIden( ff ) = strcat( lDir( ff ).name( 1:(end-4) ) );
    subjPath( ff ) = strcat( lDir( ff ).folder );

end

%% Determine whether to scale b1p
if scaleb1p

    meanb1pSubjOverall = mean( meanb1pSubj );

    for ff = 1:numFM
        b1p_CP( :, :, :, ff ) = b1p_CP( :, :, :, ff ) * ( meanb1pSubjOverall / meanb1pSubj( ff ) );
        b1p( :, :, :, :, ff ) = b1p( :, :, :, :, ff ) * ( meanb1pSubjOverall / meanb1pSubj( ff ) );
    end
end

%% Save Database
savSt = struct;
savSt.x = xn;
savSt.y = yn;
savSt.z = zn;
savSt.X = Xn;
savSt.Y = Yn;
savSt.Z = Zn;
savSt.b1p = b1p;
savSt.b1p_CP = b1p_CP;
savSt.db0 = db0;
savSt.db0shim = db0shim;
savSt.roi_body = roi_body;
savSt.roi_brain = roi_brain;
savSt.subjIden = subjIden;
savSt.subjPath = subjPath;
savSt.numSubj = numFM;
savSt.zi = zi;

%% Save database
if ~isfolder( saveDirPath )
    mkdir( saveDirPath )
end

%% Save fieldmaps
try
    save( fullfile( saveDirPath, saveName ), '-struct',...
        'savSt' );
catch
    save( fullfile( saveDirPath, saveName ), '-struct',...
        'savSt', "-v7.3" );
end

%% Change warning back
toobigWarn.state = 'on';
toobigWarn;

%% Plot and save figs
for ff = 1:numFM

    b1p_CP_iter = b1p_CP( :, :, :, ff );
    db0_iter = db0( :, :, :, ff );
    db0shim_iter = db0shim( :, :, :, ff );
    brainMask_iter = roi_brain( :, :, :, ff );
    bodyMask_iter = roi_body( :, :, :, ff );
    weightMask_iter = weightMask( :, :, :, ff );
    zi_opt_iter = zi{ ff }; 
    
    Fields = struct;
    Fields.gyro = 267.5e6; % rad / ( T sec )
    Fields.b1p = b1p_CP_iter;
    Fields.DB0 = db0_iter;
    Fields.roi_brain = brainMask_iter;
    Fields.roi_body = bodyMask_iter;

    % Get some names
    titleText = sprintf( "Subj: %s", subjIden( ff ) );
    weightFigName = sprintf( 'weight_subj_%s', subjIden( ff ) );
    b1pFigName = sprintf( 'b1p_subj_%s', subjIden( ff ) );
    db0FigName = sprintf( 'db0_subj_%s', subjIden( ff ) );

    % Plot Figs
    [ db0ArrayFig, b1pArrayFig, brainMaskArrayFig, bodyMaskArrayFig ] =...
        generateFieldMapImgs( Fields, binVis );

    % Plot Shim Maps
    pause( 0.5 );
    shimFig = plotShimFigs( db0_iter, db0shim_iter, zi_opt_iter, brainMask_iter, titleText, binVis );

    % Plot b1p Maps
    pause( 0.5 );
    b1pFig = plotb1pFigs( b1p_CP_iter, zi_opt_iter, brainMask_iter, titleText, binVis );

    % b1pCoilFigs = cell( Adj.coils, 1 );
    % b1pCoilFigNames = strings( Adj.coils, 1 );
    % 
    % for cc = 1:Adj.coils
    %     coilTitleText = sprintf( "Subj: %s, Coil: %i", subjIden( ff ), cc );
    %     b1pCoilFigs{ cc } = plotb1pFigs( squeeze( b1p_iter( cc, :, :, : ) ), zi_opt_iter, brainMask_iter, coilTitleText, binVis );
    %     b1pCoilFigNames( cc, 1 ) = sprintf( 'b1p_subj_%s_coil_%i', subjIden( ff ), cc );
    % end

    % Plot Weight Fig
    pause( 0.5 );
    weightFig = plotWeightFig( weightMask_iter, zi_opt_iter, brainMask_iter, titleText, binVis );

    %% Save Figs
    if saveFigs

        saveDirPathIter = fullfile( saveDirPath, subjIden( ff ) );
        if ~isfolder( saveDirPathIter )
            mkdir( saveDirPathIter )
        end

        if ff == 1

            fileRecordInit = getFunctionRecord( currFile );
            fileRecordSave = [ strcat( "% ", currFile ); ""; ""; fileRecordInit ];
            saveTxtFromStrArray( fullfile( saveDirPath, 'fileRecord.txt' ), fileRecordSave );

            generateMaskRecordFile = fullfile( currDir, "func", "generateMaskFromB0.m" );
            generateMaskRecordInit = getFunctionRecord( generateMaskRecordFile );
            generateMaskRecordSave = [ strcat( "% ", generateMaskRecordFile ); ""; ""; generateMaskRecordInit ];
            saveTxtFromStrArray( fullfile( saveDirPath, 'generateMaskFromB0.txt' ), generateMaskRecordSave );

        end

        saveas( db0ArrayFig, fullfile( saveDirPathIter, "db0_array.png" ) );
        saveas( db0ArrayFig, fullfile( saveDirPathIter, "db0_array.fig" ) );

        saveas( b1pArrayFig, fullfile( saveDirPathIter, "b1p_array.png" ) );
        saveas( b1pArrayFig, fullfile( saveDirPathIter, "b1p_array.fig" ) );

        saveas( brainMaskArrayFig, fullfile( saveDirPathIter, "brainMask_array.png" ) );
        saveas( brainMaskArrayFig, fullfile( saveDirPathIter, "brainMask_array.fig" ) );

        saveas( bodyMaskArrayFig, fullfile( saveDirPathIter, "bodyMask_array.png" ) );
        saveas( bodyMaskArrayFig, fullfile( saveDirPathIter, "bodyMask_array.fig" ) );

        saveas( shimFig, fullfile( saveDirPathIter, strcat( db0FigName, ".fig" ) ) );
        saveas( shimFig, fullfile( saveDirPathIter, strcat( db0FigName, ".png" ) ) );

        saveas( b1pFig, fullfile( saveDirPathIter, strcat( b1pFigName, ".fig" ) ) );
        saveas( b1pFig, fullfile( saveDirPathIter, strcat( b1pFigName, ".png" ) ) );

        % for cc = 1:Adj.coils
        %     saveas( b1pCoilFigs{ cc }, fullfile( saveDirPathIter, strcat( b1pCoilFigNames( cc ), ".fig" ) ) );
        %     saveas( b1pCoilFigs{ cc }, fullfile( saveDirPathIter, strcat( b1pCoilFigNames( cc ), ".png" ) ) );
        % end

        saveas( weightFig, fullfile( saveDirPathIter, strcat( weightFigName, ".fig" ) ) );
        saveas( weightFig, fullfile( saveDirPathIter, strcat( weightFigName, ".png" ) ) );

        %% Close figures
        close( [...
            db0ArrayFig, b1pArrayFig, brainMaskArrayFig, bodyMaskArrayFig, ...
            shimFig, b1pFig, weightFig ] );

    end
end

%% Helper functions
% ----------------------------------------------------------------------- %
function weightFig = plotWeightFig( weightArray, zi, roi_brain, titleText, binVis )

% pre-image computation
percHorzReduce = 0.20;
percVertReduce = 0.20;

numImgs = 56;
sliceIdx = sliceIdxNumImgs( numImgs, roi_brain( :, :, zi(1):zi(end) ) );

axFontSize = 22;
titleFontSize = 26;
outlineLineWidth = 2.00;
cmap_field = colorcet( 'L04' );
cAx_field = [ 0, 1 ];
contCol = 'w';

b1p_planes = rearrangeCutImgStack( weightArray(:,:, zi(1) - 1 + sliceIdx), percHorzReduce, percVertReduce );
roi_planes = rearrangeCutImgStack( roi_brain( :, :, zi(1) - 1 + sliceIdx ), percHorzReduce, percVertReduce );

figSize = [ 1 1 1180 1000 ];

if binVis
    weightFig = figure('color', 'white',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    weightFig = figure('color', 'white',...
        'Visible','off','position', figSize);
end

tl = tiledlayout( 1, 1, 'padding', 'compact', 'tilespacing', 'compact' );
pause(0.5)
ax1 = nexttile( 1 );
hold( ax1, 'on' );
imagesc( ax1, fliplr( permute( b1p_planes, [2, 1] ) ) );
set( ax1, 'ydir', 'normal' );
axis( ax1, 'image' );
xticklabels( ax1, [] );
yticklabels( ax1, [] );
clim( ax1, cAx_field );
colormap( ax1, cmap_field );
[~, contMagAx] =  imcontour( fliplr(permute( roi_planes, [2, 1] )), 1);
contMagAx.Color = contCol;
contMagAx.LineWidth = outlineLineWidth;
cbMag = colorbar(ax1);
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = axFontSize;
ylabel(cbMag, 'Weight [a.u.]', 'Interpreter','latex');
title( ax1, 'Weight Array',...
    'interpreter', 'latex', 'fontsize', titleFontSize );
pause(0.2);

title( tl,...
    strcat( "\textbf{", insertBefore(titleText, "_", "\"), " Weight Array for Brain Mask", "}" ),...
    'fontsize', titleFontSize+3, 'interpreter', 'latex' );

if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function smoothArray = smoothKernel( mArray )
sigma = 1 * ones(3,1);
filtSize = 5;
smoothArray = imgaussfilt3( mArray, sigma, 'filtersize', filtSize );
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function fieldFig = plotShimFigs( db0, db0shim, zi, roi_brain, titleText, binVis )

% pre-image computation
percHorzReduce = 0.00;
percVertReduce = 0.00;

gyro = 267.5e6;
gyrobar = gyro/(2*pi);

numImgs = 56;
sliceIdx = sliceIdxNumImgs( numImgs, roi_brain( :, :, zi(1):zi(end) ) );

axFontSize = 22;
titleFontSize = 26;
outlineLineWidth = 1.75;
cmap_field= colorcet( 'D01' );
cAx_field = [ -100, 100 ]; % Hertz
contCol = 'k';

db0_planes = gyrobar * rearrangeCutImgStack( db0(:,:, zi(1) - 1 + sliceIdx), percHorzReduce, percVertReduce );
db0shim_planes = gyrobar * rearrangeCutImgStack( db0shim(:,:, zi(1) - 1 + sliceIdx), percHorzReduce, percVertReduce );
roi_planes = rearrangeCutImgStack( roi_brain( :, :, zi(1) - 1 + sliceIdx ), percHorzReduce, percVertReduce );

figSize = [ 1 1 1270 550 ];

if binVis
    fieldFig = figure('color', 'white',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    fieldFig = figure('color', 'white',...
        'Visible','off','position', figSize);
end

tl = tiledlayout( 1, 2, 'padding', 'compact', 'tilespacing', 'compact' );

pause(0.5)
ax1 = nexttile( 1 );
hold( ax1, 'on' );
imagesc( ax1, fliplr( permute( db0_planes, [2, 1] ) ) );
set( ax1, 'ydir', 'normal' );
axis( ax1, 'image' );
xticklabels( ax1, [] );
yticklabels( ax1, [] );
clim( ax1, cAx_field );
colormap( ax1, cmap_field );
[~, contMagAx] =  imcontour( fliplr(permute( roi_planes, [2, 1] )), 1);
contMagAx.Color = contCol;
contMagAx.LineWidth = outlineLineWidth;
title( ax1, '\textbf{ Unshimmed }', 'interpreter', 'latex',...
     'fontsize', titleFontSize);
pause(0.2);

ax2 = nexttile( 2 );
hold( ax2, 'on' );
imagesc( ax2, fliplr( permute( db0shim_planes, [2, 1] ) ) );
set( ax2, 'ydir', 'normal' );
axis( ax2, 'image' );
xticklabels( ax2, [] );
yticklabels( ax2, [] );
clim( ax2, cAx_field );
colormap( ax2, cmap_field );
[~, contMagAx] =  imcontour( fliplr(permute( roi_planes, [2, 1] )), 1);
contMagAx.Color = contCol;
contMagAx.LineWidth = outlineLineWidth;
cbMag = colorbar(ax2);
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = axFontSize;
ylabel(cbMag, 'Field [Hz]', 'Interpreter','latex');
title( ax2, '\textbf{ SH2 Shim }',...
    'interpreter', 'latex', 'fontsize', titleFontSize );
pause(0.2);

title( tl,...
    strcat( "\textbf{", insertBefore(titleText, "_", "\"), "}" ),...
    'fontsize', titleFontSize+3, 'interpreter', 'latex' );

if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function fieldFig = plotb1pFigs( b1p, zi, roi_brain, titleText, binVis )

b1p = abs( b1p );

% pre-image computation
percHorzReduce = 0.00;
percVertReduce = 0.00;

numImgs = 56;
sliceIdx = sliceIdxNumImgs( numImgs, roi_brain( :, :, zi(1):zi(end) ) );

axFontSize = 22;
titleFontSize = 26;
outlineLineWidth = 1.75;
cmap_field = colorcet( 'L08' );
cAx_field = [ 0, 175 ]; % nT/V
contCol = 'w';

b1p_planes = rearrangeCutImgStack( (1e9) * b1p(:,:, zi(1) - 1 + sliceIdx), percHorzReduce, percVertReduce );
roi_planes = rearrangeCutImgStack( roi_brain( :, :, zi(1) - 1 + sliceIdx ), percHorzReduce, percVertReduce );

figSize = [ 1 1 635 550 ];

if binVis
    fieldFig = figure('color', 'white',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    fieldFig = figure('color', 'white',...
        'Visible','off','position', figSize);
end

tl = tiledlayout( 1, 1, 'padding', 'compact', 'tilespacing', 'compact' );
pause(0.5)
ax1 = nexttile( 1 );
hold( ax1, 'on' );
imagesc( ax1, fliplr( permute( b1p_planes, [2, 1] ) ) );
set( ax1, 'ydir', 'normal' );
axis( ax1, 'image' );
xticklabels( ax1, [] );
yticklabels( ax1, [] );
clim( ax1, cAx_field );
colormap( ax1, cmap_field );
[~, contMagAx] =  imcontour( fliplr(permute( roi_planes, [2, 1] )), 1);
contMagAx.Color = contCol;
contMagAx.LineWidth = outlineLineWidth;
cbMag = colorbar(ax1);
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = axFontSize;
ylabel(cbMag, 'Field [nT/V]', 'Interpreter','latex');
title( ax1, '$B_1^+$',...
    'interpreter', 'latex', 'fontsize', titleFontSize );
pause(0.2);

title( tl,...
    strcat( "\textbf{", insertBefore(titleText, "_", "\"), " CP Mode", "}" ),...
    'fontsize', titleFontSize+3, 'interpreter', 'latex' );

if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function scalarInterpVec = interpScalar( X, Y, Z, scalarArray, r_newpos_old )
extrapval = 0;
method = 'spline';
scalarInterpVec = interpn(...
    X, Y, Z, scalarArray,...
    r_newpos_old(:,1), r_newpos_old(:,2), r_newpos_old(:,3),...
    method, extrapval );
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function db0shim = calcSH2DB0Correction( db0, roi_brain, roi_body, X, Y, Z )

[ nx, ny, nz ] = size( db0 );

%% Generate 2nd Order Spherical Harmonics Shape

nSH = 2;
SH_fields = zeros( nx, ny, nz, (nSH^2+2*nSH+1) );

ntot = 0;
for nn = 0 : 2 % spherical harmonics convention
    for mm = -nn : nn % spherical harmonics convention
        ntot = ntot + 1;
        field = leg_rec_harmonic_cz( nn, mm, X, Y, Z );
        SH_fields( : , : , : , ntot ) = field;
    end
end

SH_fields = permute(SH_fields, [4, 1, 2, 3]);

%% Create Vector of Spherical Harmonics in the ROI
A_SH = ( SH_fields( : , roi_brain ) ).';
max_A_SH_per_harmonic = ( max( abs( A_SH ), [], 1 ) ).';
A_SH = A_SH ./ max_A_SH_per_harmonic.';
SH_fields_scaling = permute( repmat( reshape( max_A_SH_per_harmonic, [ 1, 1, 1, length(max_A_SH_per_harmonic ) ] ), size( roi_brain ) ), [ 4, 1, 2, 3 ] );
SH_fields = SH_fields ./ SH_fields_scaling;

%% Generate target field for shimming
db0vec = db0( roi_brain );

%% solve least-squares problem
shimCurr = lsqminnorm( A_SH, -db0vec );

%% create new fieldmap
db0shim = db0 + squeeze( sum( shimCurr .* SH_fields, 1 ) );
db0shim = db0shim .* double( roi_body );

end
% ----------------------------------------------------------------------- %

