function [ spokesTiming ] = determineSpokesDuration( pulse, gradMax, gradSlewRate, roundTime )

if nargin < 4
    % round to nearest time
    roundTime = 5e-6;
end

numDigits = 13;

%% Initialize timing struct and get num spokes
spokesTiming = struct;

numSpokes = pulse.numSpokes;
spokesTiming.numSpokes = numSpokes;

if isfield( pulse, 'gyro' )
    spokesTiming.gyro = pulse.gyro;
else
    pulse.gyro = 267.5e6;
    spokesTiming.gyro = pulse.gyro;
end

%% Get ramp time
rampTimeMaxGrad = roundUpToNearestMicrosecond( gradMax / gradSlewRate, roundTime, numDigits );
% kSpaceTraversalDuringMaxGradRamp = 0.5 * rampTimeMaxGrad * gradMax * pulse.gyro/ ( 2 * pi );

%% Determine lengths of the central spoke pulses
% First start with the last spoke, which we will assume is at kx = 0, ky = 0

if isfield( pulse, 'centralSpokeTBW' )
    centralSpokeTBW = pulse.centralSpokeTBW;
elseif isfield( pulse, 'spokeTBW' )
    centralSpokeTBW = pulse.spokeTBW;
else
    error( "Unknown central spoke time-bandwidth product" );
end
if isfield( pulse, 'centralSpokeLength' )
    centralSpokeLength = roundToNearestMicrosecond( pulse.centralSpokeLength, roundTime, numDigits );
elseif isfield( pulse, 'spokeLength' )
    centralSpokeLength = roundToNearestMicrosecond( pulse.spokeLength, roundTime, numDigits );
else
    error( "Unknown central spoke length." )
end

% Get the gradient time product pg. 117 Nishimura. End of section 6.2
GssTauCentralSpoke = 2 * pi * centralSpokeTBW / ( pulse.gyro * pulse.sliceThickness );

% Determine if we can achieve this gradient time product with the current
% pulse length
constTimeMaxGradCentralSpoke = roundUpToNearestMicrosecond( max( centralSpokeLength - 2*rampTimeMaxGrad, 0 ), roundTime, numDigits );

if ( round( (0.5*2*rampTimeMaxGrad*gradMax + constTimeMaxGradCentralSpoke * gradMax), numDigits ) < round( GssTauCentralSpoke, numDigits ) )
        
    % can't achieve gradient integral
    centralSpokeSlewTime = rampTimeMaxGrad;
    centralSpokeConstantTime = roundUpToNearestMicrosecond(...
        (GssTauCentralSpoke - 2*rampTimeMaxGrad*gradMax*0.5)/gradMax, roundTime, numDigits );
    centralSpokeLength = centralSpokeConstantTime + 2 * centralSpokeSlewTime;
    centralSpokeGradMag = GssTauCentralSpoke / ( centralSpokeConstantTime + centralSpokeSlewTime );

    warning( "Need to elongate central spoke timing" );
    warning( "centralSpokeLength is now %i us", round( centralSpokeLength * 1e6 ) );
elseif ( round( gradSlewRate * ( 0.5 * centralSpokeLength )^2, numDigits ) < round( GssTauCentralSpoke, numDigits ) )
    % can't achieve gradient integral, but can during slew
    centralSpokeLength = roundUpToNearestMicrosecond( ...
        sqrt( GssTauCentralSpoke / gradSlewRate ) / 0.5, 2 * roundTime, numDigits );
    centralSpokeSlewTime = centralSpokeLength / 2;
    centralSpokeConstantTime = 0;
    centralSpokeGradMag = GssTauCentralSpoke / centralSpokeSlewTime;

    warning( "Need to elongate central spoke timing" );
    warning( "centralSpokeLength is now %i us", round( centralSpokeLength * 1e6 ) );
else 
    % can achieve gradient integral. let's find out lengths and values
    possibleGradCentralSpokeValues = roots( [ -1/gradSlewRate, centralSpokeLength, -GssTauCentralSpoke ] );
    if isreal( possibleGradCentralSpokeValues )
        centralSpokeGradMag = possibleGradCentralSpokeValues(...
            ( possibleGradCentralSpokeValues > 0 ) &...
            ( possibleGradCentralSpokeValues < gradMax ) );
        centralSpokeSlewTime = roundUpToNearestMicrosecond( centralSpokeGradMag/gradSlewRate, roundTime, numDigits );
        centralSpokeConstantTime = roundUpToNearestMicrosecond(...
            centralSpokeLength - 2 * centralSpokeSlewTime, roundTime, numDigits );
        centralSpokeLength = centralSpokeConstantTime + 2 * centralSpokeSlewTime;
        % re-adjust central spoke grad value based on the new timing
        centralSpokeGradMag = GssTauCentralSpoke / ( centralSpokeConstantTime + centralSpokeSlewTime );
    else
        error( "Returned complex central spoke grad values." )
    end

end

spokesTiming.centralSpokeSlewTime = centralSpokeSlewTime;
spokesTiming.centralSpokeConstantTime = centralSpokeConstantTime;
spokesTiming.centralSpokeLength = centralSpokeLength;
spokesTiming.centralSpokeGradMag = centralSpokeGradMag;
spokesTiming.centralSpokeTBW = centralSpokeTBW;

%% determine rephase time
if round( ( ( 0.5 * GssTauCentralSpoke ) / gradMax ), numDigits ) < rampTimeMaxGrad
    possibleRephaseTimingValues = roots( [ -gradSlewRate, 0, 0.5 * GssTauCentralSpoke ] );
    possibleRephaseTimingValues = ...
        possibleRephaseTimingValues( possibleRephaseTimingValues > 0 ); 
    rephaseSlewTime = roundUpToNearestMicrosecond(min( possibleRephaseTimingValues ), roundTime, numDigits );
    rephaseConstantTime = 0;
    rephaseLength = 2 * rephaseSlewTime + rephaseConstantTime;
    rephaseGradMag = 0.5 * GssTauCentralSpoke / rephaseSlewTime;
else
    rephaseSlewTime = rampTimeMaxGrad;
    rephaseConstantTime = roundUpToNearestMicrosecond(...
        ( 0.5 * GssTauCentralSpoke - 2*rampTimeMaxGrad*gradMax*0.5)/gradMax, roundTime, numDigits );
    rephaseLength = 2 * rephaseSlewTime + rephaseConstantTime;
    rephaseGradMag = 0.5 * GssTauCentralSpoke / ( rephaseConstantTime + rephaseSlewTime );
end

spokesTiming.rephaseSlewTime = rephaseSlewTime;
spokesTiming.rephaseConstantTime = rephaseConstantTime;
spokesTiming.rephaseLength = rephaseLength;
spokesTiming.rephaseGradMag = rephaseGradMag;

%% Determine lengths of the non-central spoke pulses
if numSpokes > 1

    nonCentralSpokeGradMag = centralSpokeGradMag;

    nonCentralSpokeTBWSpecified = isfield( pulse, 'nonCentralSpokeTBW' ) | isfield( pulse, 'spokeTBW' ) ;
    nonCentralSpokeLengthSpecified = isfield( pulse, 'nonCentralSpokeLength' ) | isfield( pulse, 'spokeLength' ) ;

    if ( nonCentralSpokeTBWSpecified ) && ( nonCentralSpokeLengthSpecified )
        warning( "nonCentralSpokeTBW and nonCentralSpokeLength both specified. Default to TBW." )
        nonCentralSpokeLengthSpecified = false;
    end

    if nonCentralSpokeTBWSpecified

        if isfield( pulse, 'nonCentralSpokeTBW' )
            nonCentralSpokeTBW = pulse.nonCentralSpokeTBW;
        elseif isfield( pulse, 'spokeTBW' )
            nonCentralSpokeTBW = pulse.spokeTBW;
        end

        % given the TBW, we know that the gradient has to be the same as
        % for the central slice to keep the slice location in the same
        % position, so we can figure out the timing

        if round( ( 2 * 0.5 * centralSpokeSlewTime * nonCentralSpokeGradMag ), numDigits ) >...
                round( ( nonCentralSpokeTBW / ( pulse.gyro / ( 2*pi ) * pulse.sliceThickness ) ), numDigits )
            
            % Need to increase TBW
            nonCentralSpokeTBW = 2 * 0.5 * centralSpokeSlewTime * nonCentralSpokeGradMag * pulse.sliceThickness * pulse.gyro / ( 2*pi );
            nonCentralSpokeSlewTime = centralSpokeSlewTime;
            nonCentralSpokeConstantTime = 0;
            nonCentralSpokeLength = 2 * nonCentralSpokeSlewTime;

            warning( "Need to increase non-central spoke TBW" );
            warning( "nonCentralSpokeTBW is now %g", nonCentralSpokeTBW );
            
        else
            nonCentralSpokeSlewTime = centralSpokeSlewTime;
            nonCentralSpokeConstantTime = nonCentralSpokeTBW / ( pulse.gyro / ( 2*pi ) * pulse.sliceThickness * nonCentralSpokeGradMag )...
                - nonCentralSpokeSlewTime;
            nonCentralSpokeConstantTime = roundUpToNearestMicrosecond( nonCentralSpokeConstantTime, roundTime, numDigits );
            nonCentralSpokeLength = nonCentralSpokeConstantTime + 2 * nonCentralSpokeSlewTime;
            nonCentralSpokeTBW = ( nonCentralSpokeConstantTime + nonCentralSpokeSlewTime ) * nonCentralSpokeGradMag * pulse.sliceThickness * pulse.gyro / ( 2*pi );
        end

    elseif nonCentralSpokeLengthSpecified

        if isfield( pulse, 'nonCentralSpokeLength' )
            nonCentralSpokeLength = roundToNearestMicrosecond( pulse.nonCentralSpokeLength, roundTime, numDigits );
        elseif isfield( pulse, 'spokeLength' )
            nonCentralSpokeLength = roundToNearestMicrosecond( pulse.spokeLength, roundTime, numDigits );
        end

        if ( 2*centralSpokeSlewTime ) > nonCentralSpokeLength
            nonCentralSpokeLength = 2*centralSpokeSlewTime;
            nonCentralSpokeSlewTime = centralSpokeSlewTime;
            nonCentralSpokeConstantTime = 0;

            nonCentralSpokeTBW = pulse.gyro / ( 2*pi ) * centralSpokeGradMag * ( nonCentralSpokeSlewTime ) * pulse.sliceThickness;

            warning( "Need to elongate non-central spoke timing" );
            warning( "nonCentralSpokeLength is now %i us", round( nonCentralSpokeLength * 1e6 ) );
        else
            nonCentralSpokeSlewTime = centralSpokeSlewTime;
            nonCentralSpokeConstantTime = nonCentralSpokeLength - 2*nonCentralSpokeSlewTime;
            nonCentralSpokeTBW = pulse.gyro / ( 2*pi ) * centralSpokeGradMag * ( nonCentralSpokeConstantTime + nonCentralSpokeSlewTime ) * pulse.sliceThickness;
        end
        
    else
        error( "Unknown non-central spoke time-bandwidth product or spoke length" );
    end

else
    nonCentralSpokeLength = 0;
    nonCentralSpokeSlewTime = 0;
    nonCentralSpokeConstantTime = 0;
    nonCentralSpokeGradMag = 0;
    nonCentralSpokeTBW = 0;
end

spokesTiming.nonCentralSpokeSlewTime = nonCentralSpokeSlewTime;
spokesTiming.nonCentralSpokeConstantTime = nonCentralSpokeConstantTime;
spokesTiming.nonCentralSpokeLength = nonCentralSpokeLength;
spokesTiming.nonCentralSpokeGradMag = nonCentralSpokeGradMag;
spokesTiming.nonCentralSpokeTBW = nonCentralSpokeTBW;

%% Set gradient blip length between spokes
if numSpokes > 1
    gradMaxBlip = min( gradMax, 10e-3 );
    spokeBlipLength = 2 * roundUpToNearestMicrosecond( gradMaxBlip / gradSlewRate, roundTime, numDigits );

    spokesTiming.spokeBlipLength = spokeBlipLength;
else
    spokesTiming.spokeBlipLength = 0;
end

%% Determine necessary Grad slice select slew before the central spoke
% need to ensure that the non-central spokes are centered about zero
if centralSpokeTBW ~= nonCentralSpokeTBW
    kmax_centralSpoke = pulse.gyro/(2*pi)*( rephaseSlewTime + rephaseConstantTime ) * rephaseGradMag;
    kmax_nonCentralSpoke = 0.5 * pulse.gyro/(2*pi)*( nonCentralSpokeSlewTime + nonCentralSpokeConstantTime ) * nonCentralSpokeGradMag;

    finalBlipGradTimeArea = ( kmax_nonCentralSpoke - kmax_centralSpoke ) / ( -pulse.gyro / ( 2 * pi ) );
    absFinalBlipGradTimeArea = abs( finalBlipGradTimeArea );

    if round( ( absFinalBlipGradTimeArea / gradMax ), numDigits ) < rampTimeMaxGrad
        possibleLastBlipTimingValues = roots( [ -gradSlewRate, 0, absFinalBlipGradTimeArea ] );
        if isreal( possibleLastBlipTimingValues )
            possibleLastBlipTimingValues = ...
                possibleLastBlipTimingValues( possibleLastBlipTimingValues > 0 );
            lastBlipSlewTime = roundUpToNearestMicrosecond(min( possibleLastBlipTimingValues ), roundTime, numDigits );
            lastBlipConstantTime = 0;
            lastBlipLength = 2 * lastBlipSlewTime + lastBlipConstantTime;
            lastBlipGradMag = absFinalBlipGradTimeArea / lastBlipSlewTime;
        else
            error( "Returned negative last blip timing values." )
        end
    else
        lastBlipSlewTime = rampTimeMaxGrad;
        lastBlipConstantTime = roundUpToNearestMicrosecond(...
            ( absFinalBlipGradTimeArea - 2*rampTimeMaxGrad*gradMax*0.5)/gradMax, roundTime, numDigits );
        lastBlipLength = 2 * lastBlipSlewTime + lastBlipConstantTime;
        lastBlipGradMag = absFinalBlipGradTimeArea / ( lastBlipConstantTime + lastBlipSlewTime );
    end

    if finalBlipGradTimeArea > 0
        lastBlipGradMag = -1 * lastBlipGradMag;
    end

else
    if numSpokes == 1
        lastBlipLength = 0;
        lastBlipGradMag = 0;
        lastBlipSlewTime = 0;
        lastBlipConstantTime = 0;
    else
        lastBlipLength = spokeBlipLength;
        lastBlipGradMag = 0;
        lastBlipSlewTime = spokeBlipLength / 2;
        lastBlipConstantTime = 0;

    end
end

spokesTiming.lastBlipSlewTime = lastBlipSlewTime;
spokesTiming.lastBlipConstantTime = lastBlipConstantTime;
spokesTiming.lastBlipLength = lastBlipLength;
spokesTiming.lastBlipGradMag = lastBlipGradMag;

spokesTiming.numSlices = 1;

%% Determine amplitude modulation frequency
spokesTiming.modFreq = pulse.gyro/(2*pi) * spokesTiming.centralSpokeGradMag * pulse.sliceLocation; % in Hz

%% Determine total pulse length
if numSpokes > 1
    pulseLength = rephaseLength + centralSpokeLength + lastBlipLength +...
        ( numSpokes - 2 ) * spokeBlipLength + ( numSpokes - 1 ) * nonCentralSpokeLength;
else
    pulseLength = rephaseLength + centralSpokeLength;
end

spokesTiming.pulseLength = pulseLength;

%% Get slice information
if ~isfield( pulse, 'sliceDirection' )
    pulse.sliceDirection = [ 0; 0; 1 ];
end
spokesTiming.sliceDirection = pulse.sliceDirection;
spokesTiming.sliceLocation = pulse.sliceLocation;
spokesTiming.sliceThickness = pulse.sliceThickness;

spokesTiming.sliceBounds = determineSliceBounds( spokesTiming.sliceLocation, spokesTiming.sliceThickness );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function out = roundToNearestMicrosecond( num, roundTime, numDigits )
if nargin < 3
    numDigits = 13;
end
out = roundTime * floor( round( num/roundTime, numDigits ) );
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function out = roundUpToNearestMicrosecond( num, roundTime, numDigits )
if nargin < 3
    numDigits = 13;
end
out = roundTime * ceil( round( num/roundTime, numDigits ) );
end
% ----------------------------------------------------------------------- %