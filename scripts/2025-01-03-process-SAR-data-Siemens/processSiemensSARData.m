%% File Initialization
restoredefaultpath;
clearvars;
% close all;
home;

%% Description
% This script will load in the SAR provided for the Siemens Terra.X
% and into a database that will be used for pulse optimization
%
% Script was prepared to process the following files:
% 
% SarDataUserTerraX.mat
% SysDataUserTerraX.mat

%% Set Paths
currFile = strcat( mfilename( 'fullpath' ), ".m" );
currDir = fileparts( currFile );

cld = strsplit( currDir, 'scripts' );
cld = cld{ 1 };
addpath( genpath( fullfile( cld, 'util' ) ) );
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

%% Set save dir
saveDirPath = fullfile( tld, 'data', 'fields', '2025-01-02-Siemens-Nova-8ch-pTx-SAR' );
saveName = "Nova_8ch_pTx_SAR.mat";

%% Change warning temporarily
toobigWarn = warning('error', 'MATLAB:save:sizeTooBigForMATFile');

%% List Directory for loading
sysDataFile = fullfile( saveDirPath, "SarSysData", "SysDataUserTerraX.mat" );
sarDataFile = fullfile( saveDirPath, "SarSysData", "SarDataUserTerraX.mat" );

sys = load( sysDataFile );
sar = load( sarDataFile );

%% Process SAR matrices
gsarid = 0; % ID for global SAR matrix
pbsarid = 1; % ID for partial body SAR matrix
lsarid = 5; % ID for local SAR matrices
rfcoilenergyid = 9; % ID for RF coil
rfpaenergyid = 10; % ID for RFPA

QGlobal = sar.ZZ( :, :, sar.ZZtype == gsarid );
QGlobalLim = sar.VOPTypeLimit( sar.VOPTypeLimitID == gsarid );

QPartialBody = sar.ZZ( :, :, sar.ZZtype == pbsarid );
QPartialBodyLim = sar.VOPTypeLimit( sar.VOPTypeLimitID == pbsarid );

VOPs = sar.ZZ( :, :, sar.ZZtype == lsarid );
VOPsLim = sar.VOPTypeLimit( sar.VOPTypeLimitID == lsarid );

RFCoilEnergy = sar.ZZ( :, :, sar.ZZtype == rfcoilenergyid );
RFCoilEnergyLim = sar.VOPTypeLimit( sar.VOPTypeLimitID == rfcoilenergyid );

RFPAEnergy = sar.ZZ( :, :, sar.ZZtype == rfpaenergyid );
RFPAEnergyLim = sar.VOPTypeLimit( sar.VOPTypeLimitID == rfpaenergyid );

%% Process system limits
gradmax = sys.gmax;
gradslewmax = sys.smax;
rfmax = sys.umax;

%% Save Database
savSt = struct;
savSt.QGlobal = QGlobal;
savSt.QGlobalLim = QGlobalLim;
savSt.QPartialBody = QPartialBody;
savSt.QPartialBodyLim = QPartialBodyLim;
savSt.VOPs = VOPs;
savSt.VOPsLim = VOPsLim;
savSt.RFCoilEnergy = RFCoilEnergy;
savSt.RFCoilEnergyLim = RFCoilEnergyLim;
savSt.RFPAEnergy = RFPAEnergy;
savSt.RFPAEnergyLim = RFPAEnergyLim;

savSt.gradmax = gradmax;
savSt.gradslewmax = gradslewmax;
savSt.rfmax = rfmax;

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

%% Helper functions
% ----------------------------------------------------------------------- %


