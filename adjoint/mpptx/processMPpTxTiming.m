function [ oc, pulse, opt ] = processMPpTxTiming( oc, pulse, opt )

%% Add multiphoton parameters to opt struct
opt.fz_mpptx = pulse.fz_mpptx;
opt.wz_mpptx = 2 * pi * opt.fz_mpptx;

opt.dfxy_mpptx = pulse.dfxy_mpptx;
opt.dwxy_mpptx = 2 * pi * opt.dfxy_mpptx;

opt.gyro = 267.5e6;

%% Check if constraints are present. If not, add in the needed ones
constraintList = keys( pulse.constraints );

if ~any( contains( constraintList, "RF-max", 'ignorecase', true ) )
    pulse.constraints( "RF-max" ) = 550; % V
end
if ~any( contains( constraintList, "grad-max", 'ignorecase', true ) )
    pulse.constraints( "grad-max" ) = 25e-3; % T/m
end
if ~any( contains( constraintList, "grad-slew-rate", 'ignorecase', true ) )
    pulse.constraints( "grad-slew-rate" ) = 200; % T/m/sec
end
if ( opt.numZCoils > 0 ) && pulse.shimarray
    if ~any( contains( constraintList, "shim-max", 'ignorecase', true ) )
        pulse.constraints( "shim-max" ) = 50; % Amp-turns
    end
    if ~any( contains( constraintList, "shim-slew-rate", 'ignorecase', true ) )
        pulse.constraints( "shim-slew-rate" ) = 3e6; % Amp-turns/sec
    end
end

%% See if blip constraints are present, and if not, add blip constraints 
if ~any( contains( constraintList, "grad-blip-max", 'ignorecase', true ) )
    pulse.constraints( "grad-blip-max" ) = pulse.constraints( "grad-max" );
end
if ( opt.numZCoils > 0 ) && pulse.shimarray
    if ~any( contains( constraintList, "shim-blip-max", 'ignorecase', true ) )
        pulse.constraints( "shim-blip-max" ) = pulse.constraints( "shim-max" );
    end
end

clear constraintList;

%% Determine slew limits
numSigDigitsRound = 6;

dt = opt.dt;
tvec = opt.tvec;
pulseLength = opt.pulseLength;

dt_tol = 1e-6;

% need to determine the slew time for the RF during the RF periods of kTP
if isfield( pulse, "minRFSlewTime" )
    minRFSlewTime = pulse.minRFSlewTime;
else
    minRFSlewTime = 10e-6; % 10 us minimum slew time
end
if ~isfield( pulse, 'RFSlewTime' )
    RFSlewTime = max( [ abs( pulse.constraints( 'RF-max' ) / pulse.constraints( 'RF-slew-rate' ) ), minRFSlewTime ] );
    RFSlewTime = round( round( RFSlewTime/dt ) * dt, numSigDigitsRound, 'significant');
else
    RFSlewTime = pulse.minRFSlewTime;
end

if isfield( pulse, "minSlewtime" )
    minSlewTime = pulse.minSlewTime;
else
    minSlewTime = 10e-6; % 10 us minimum slew time
end

% get constraints
constraintList = keys( pulse.constraints );

% first need to determine maximum magnitude during the MPSP due to the slew
% rate.
if ( opt.numZCoils > 0 ) && pulse.shimarray
    maxShimMagfromMPSPSlewRate = pulse.constraints( 'shim-slew-rate' ) / opt.wz_mpptx;
    if any( contains( constraintList, "shim-max-MPSP", "ignorecase", true ) )
        shimMaxMPSP = pulse.constraints( 'shim-max-MPSP' );
    elseif any( contains( constraintList, "shim-max", "ignorecase", true ) )
        shimMaxMPSP = pulse.constraints( 'shim-max' );
    end
    pulse.constraints( 'shim-max-MPSP' ) = ...
        min( [ maxShimMagfromMPSPSlewRate, shimMaxMPSP ] );
end

maxGradMagfromMPSPSlewRate = pulse.constraints( 'grad-slew-rate' ) / opt.wz_mpptx;
if any( contains( constraintList, "grad-max-MPSP", "ignorecase", true ) )
    gradMaxMPSP = pulse.constraints( 'grad-max-MPSP' );
elseif any( contains( constraintList, "grad-max", "ignorecase", true ) )
    gradMaxMPSP = pulse.constraints( 'grad-max' );
end

pulse.constraints( 'grad-max-MPSP' ) = ...
    min( [ maxGradMagfromMPSPSlewRate, gradMaxMPSP ] );

% Determine slew times
if ( opt.numZCoils > 0 ) && pulse.shimarray
    if ~isfield( pulse, 'BlipShimSlewTime' )
        BlipShimSlewTime = max( [...
            abs( pulse.constraints( 'shim-blip-max' ) / pulse.constraints( 'shim-slew-rate' ) ),...
            minSlewTime ] );
        BlipShimSlewTime = round( round( BlipShimSlewTime/dt ) * dt, numSigDigitsRound, 'significant' );
    else
        BlipShimSlewTime = pulse.BlipShimSlewTime;
    end

    if ~isfield( pulse, 'MPSPShimSlewTime' )
        MPSPShimSlewTime = max( [ abs( pulse.constraints( 'shim-max-MPSP' ) / pulse.constraints( 'shim-slew-rate' ) ), minSlewTime ] );
        MPSPShimSlewTime = ceil( MPSPShimSlewTime/dt - 1e-6 ) * dt;
    else
        MPSPShimSlewTime = pulse.MPSPShimSlewTime;
    end
end

if ~isfield( pulse, 'BlipGradSlewTime' )
    BlipGradSlewTime = max( [ abs( pulse.constraints( 'grad-blip-max' ) / pulse.constraints( 'grad-slew-rate' ) ), minSlewTime ] );
    BlipGradSlewTime = round( round( BlipGradSlewTime/dt ) * dt, numSigDigitsRound, 'significant');
else
    BlipGradSlewTime = pulse.BlipGradSlewTime;
end

if ~isfield( pulse, 'MPSPGradSlewTime' )
    MPSPGradSlewTime = max( [ abs( pulse.constraints( 'grad-max-MPSP' ) / pulse.constraints( 'grad-slew-rate' ) ), minSlewTime ] );
    MPSPGradSlewTime = ceil( MPSPGradSlewTime/dt - 1e-6 ) * dt;
else
    MPSPGradSlewTime = pulse.MPSPGradSlewTime;
end

%% Determine points of different periods

pulse.tORSP = round( pulse.tORSP, numSigDigitsRound );
pulse.tBlip = round( pulse.tBlip, numSigDigitsRound );
pulse.tMPSP = round( pulse.tMPSP, numSigDigitsRound );

pulse.tORSP = max( [ pulse.tORSP, 2*RFSlewTime ] );
pulse.tMPSP = max( [ pulse.tMPSP, 2*RFSlewTime ] );

tend = pulse.tORSP + pulse.tBlip + pulse.tMPSP;
if abs( (tend - pulseLength)/pulseLength ) > dt_tol
    error( "Pulse length and tORSP + tBlip + tMPSP are not equivalent" );
end

tStORSP = 0;
tEndORSP = pulse.tORSP;

tStBlip = tEndORSP;
tEndBlip = tEndORSP + pulse.tBlip;

tStMPSP = tEndBlip;
tEndMPSP = tEndBlip + pulse.tMPSP;

idxtol = dt/10;
ORSP_i = find( tvec >= ( tStORSP - idxtol  ), 1, 'first' );
ORSP_f = find( tvec <= ( tEndORSP + idxtol ), 1, 'last' );
ORSP_RF_Slew_i = find( tvec >= ( ( tStORSP + RFSlewTime ) - idxtol  ), 1, 'first' );
ORSP_RF_Slew_f = find( tvec <= ( ( tEndORSP - RFSlewTime ) + idxtol ), 1, 'last' );

num_ORSP = ORSP_f - ORSP_i + 1;

Blip_i = find( tvec >= ( tStBlip - idxtol  ), 1, 'first' );
Blip_f = find( tvec <= ( tEndBlip + idxtol ), 1, 'last' );
if ( opt.numZCoils > 0 ) && pulse.shimarray
    Blip_Shim_Slew_i = find( tvec >= ( ( tStBlip + BlipShimSlewTime ) - idxtol  ), 1, 'first' );
    Blip_Shim_Slew_f = find( tvec <= ( ( tEndBlip - BlipShimSlewTime ) + idxtol ), 1, 'last' );
end

Blip_Grad_Slew_i = find( tvec >= ( ( tStBlip + BlipGradSlewTime ) - idxtol  ), 1, 'first' );
Blip_Grad_Slew_f = find( tvec <= ( ( tEndBlip - BlipGradSlewTime ) + idxtol ), 1, 'last' );

MPSP_i = find( tvec >= ( tStMPSP - idxtol  ), 1, 'first' );
MPSP_f = find( tvec <= ( tEndMPSP + idxtol ), 1, 'last' );
MPSP_RF_Slew_i = find( tvec >= ( ( tStMPSP + RFSlewTime ) - idxtol  ), 1, 'first' );
MPSP_RF_Slew_f = find( tvec <= ( ( tEndMPSP - RFSlewTime ) + idxtol ), 1, 'last' );

if ( opt.numZCoils > 0 ) && pulse.shimarray
    MPSP_Shim_Slew_i = find( tvec >= ( ( tStMPSP + MPSPShimSlewTime ) - idxtol  ), 1, 'first' );
    MPSP_Shim_Slew_f = find( tvec <= ( ( tEndMPSP - MPSPShimSlewTime) + idxtol ), 1, 'last' );
end

MPSP_Grad_Slew_i = find( tvec >= ( ( tStMPSP + MPSPGradSlewTime ) - idxtol  ), 1, 'first' );
MPSP_Grad_Slew_f = find( tvec <= ( ( tEndMPSP - MPSPGradSlewTime) + idxtol ), 1, 'last' );
num_MPSP = MPSP_f - MPSP_i + 1;

%% Exception Handling
if Blip_i == Blip_f
    Blip_i = [];
    Blip_f = [];

    if ORSP_f >= MPSP_i
        ORSP_f = ORSP_f - 1;
    end
else
    if ORSP_f >= Blip_i
        Blip_i = Blip_i + 1;
    end
    if Blip_f >= MPSP_i
        Blip_f = Blip_f - 1;
    end
end

%% Assign to opt struct
opt.ORSP_i = ORSP_i;
opt.ORSP_f = ORSP_f;
opt.ORSP_RF_Slew_i = ORSP_RF_Slew_i;
opt.ORSP_RF_Slew_f = ORSP_RF_Slew_f;
opt.num_ORSP = num_ORSP;

opt.Blip_i = Blip_i;
opt.Blip_f = Blip_f;
if ( opt.numZCoils > 0 ) && pulse.shimarray
    opt.Blip_Shim_Slew_i = Blip_Shim_Slew_i;
    opt.Blip_Shim_Slew_f = Blip_Shim_Slew_f;
end
opt.Blip_Grad_Slew_i = Blip_Grad_Slew_i;
opt.Blip_Grad_Slew_f = Blip_Grad_Slew_f;

opt.MPSP_i = MPSP_i;
opt.MPSP_f = MPSP_f;
opt.MPSP_RF_Slew_i = MPSP_RF_Slew_i;
opt.MPSP_RF_Slew_f = MPSP_RF_Slew_f;
opt.num_MPSP = num_MPSP;

if ( opt.numZCoils > 0 ) && pulse.shimarray
    opt.MPSP_Shim_Slew_i = MPSP_Shim_Slew_i;
    opt.MPSP_Shim_Slew_f = MPSP_Shim_Slew_f;
end

opt.MPSP_Grad_Slew_i = MPSP_Grad_Slew_i;
opt.MPSP_Grad_Slew_f = MPSP_Grad_Slew_f;

opt.tvec = tvec;

opt.tORSP = pulse.tORSP;
opt.tBlip = pulse.tBlip;
opt.tMPSP = pulse.tMPSP;

opt.tStORSP = tStORSP;
opt.tEndORSP = tEndORSP;
opt.tStBlip = tStBlip;
opt.tEndBlip = tEndBlip;
opt.tStMPSP = tStMPSP;
opt.tEndMPSP = tEndMPSP;

if ( opt.numZCoils > 0 ) && pulse.shimarray
    opt.BlipShimSlewTime = BlipShimSlewTime;
    opt.MPSPShimSlewTime = MPSPShimSlewTime;
end

opt.BlipGradSlewTime = BlipGradSlewTime;
opt.MPSPGradSlewTime = MPSPGradSlewTime;

if ( opt.numZCoils > 0 ) && pulse.shimarray
    pulse.BlipShimSlewTime = BlipShimSlewTime;
    pulse.MPSPShimSlewTime = MPSPShimSlewTime;
end
pulse.BlipGradSlewTime = BlipGradSlewTime;
pulse.MPSPGradSlewTime = MPSPGradSlewTime;

opt.minRFSlewTime = minRFSlewTime;
opt.minSlewTime = minSlewTime;

opt.RFSlewTime = RFSlewTime;
opt.RFSlewIntPtNum = round( RFSlewTime / dt );
opt.tORSPIntPtNum = round( opt.tORSP / dt ) - 2 * opt.RFSlewIntPtNum;
opt.tMPSPIntPtNum = round( opt.tMPSP / dt ) - 2 * opt.RFSlewIntPtNum;


end