function [ inSliceRipple, outSliceRipple, inSliceRipplePercent, outSliceRipplePercent,...
    transitionWidth, FWHM, phaseChangeAcrossSlice, magNRMSE_FA, magNRMSE_compFA,...
    xySamplePtMxVals, xySamplePtMyVals, xySamplePtFAVals,...
    xySamplePtPhVals, xySamplePtzVals, xySamplePtFAtargVals ] = ...
    calculateThruSliceMetrics( MThru, Mtarg, FAthru, FAtarg, val )

sliceLocation = val.sliceLocation;
sliceBounds = val.sliceBounds;

di_tol = abs( val.z( 2 ) - val.z( 1 ) ) / 10;
FApercentage = 1e-4;
TBWLow = 0.075;
TBWHigh = 0.925;
FWHMVal = 0.5;
% TBWFind = 0.20;
interpMethod = 'pchip';
extrap = 0;
sliceBoundExtendAmount = 1;

numSlices = length( sliceLocation );

inSliceRipple = zeros( numSlices, 1 );
outSliceRipple = zeros( numSlices, 1 );
inSliceRipplePercent = zeros( numSlices, 1 );
outSliceRipplePercent = zeros( numSlices, 1 );
transitionWidth = zeros( numSlices, 1 );
FWHM = zeros( numSlices, 1 );
phaseChangeAcrossSlice = zeros( numSlices, 1 );
xySamplePtMxVals = cell( numSlices, 1 );
xySamplePtMyVals = cell( numSlices, 1 );
xySamplePtFAtargVals = cell( numSlices, 1 );
xySamplePtFAVals = cell( numSlices, 1 );
xySamplePtPhVals = cell( numSlices, 1 );
xySamplePtzVals = cell( numSlices, 1 );

% calculate magNRMSE_complexFA
numPosSubj = size( MThru, 1 );

Mx = MThru( :, 1 );
My = MThru( :, 2 );
Mz = MThru( :, 3 );

Mzbiggerthan1 = Mz > 1;
Mzlessthanm1 = Mz < -1;

Mz( Mzbiggerthan1 ) = 1;
Mz( Mzlessthanm1 ) = -1;

subjSlicesGloIdx = val.subjSlicesGloIdxs( :, 1 );

subjSlicesLocIdx = val.subjSlicesLocIdxs{ 1, 1 };
subjNSlicesLocIdx = val.subjNSlicesLocIdxs{ 1, 1 };
numInSlice = val.numInSliceSubjs( 1, 1 );

uniqueLocIdxSlice = zeros( numInSlice, val.numUniqueXYPos );
uniqueLocIdxSliceLinIdx = sub2ind( [ numInSlice, val.numUniqueXYPos ],...
    transpose( uint32(1:numInSlice) ), val.uniqueXYPosToXYPos( subjSlicesGloIdx ) );
uniqueLocIdxSlice( uniqueLocIdxSliceLinIdx ) = ones( numInSlice, 1 );
uniqueIdxToPos = val.uniqueXYPosToXYPos( subjSlicesGloIdx );

Mx_inslice = Mx( subjSlicesLocIdx );
My_inslice = My( subjSlicesLocIdx );
Mz_inslice = Mz( subjSlicesLocIdx );
% Mx_outslice = Mx( subjNSlicesLocIdx );
% My_outslice = My( subjNSlicesLocIdx );
Mz_outslice = Mz( subjNSlicesLocIdx );

pulseFA_inslice = real( acosd( Mz_inslice ) );
pulseFA_outslice = real( acosd( Mz_outslice ) );
theta_outslice = pulseFA_outslice;

Mxyph_inslice = atan2( My_inslice, Mx_inslice );
expMxyph_inslice = exp( 1j * Mxyph_inslice );

theta_inslice = pulseFA_inslice .* expMxyph_inslice;

Mztarg = Mtarg( :, 3 );

targFA = real( acosd( Mztarg ) );
targFA_inslice = targFA( subjSlicesLocIdx );
targFA_outslice = targFA( subjNSlicesLocIdx );
thetatarg_outslice = targFA_outslice;

Mxsum_inslice_xy = transpose( sum( Mx_inslice .* uniqueLocIdxSlice, 1 ) );
Mysum_inslice_xy = transpose( sum( My_inslice .* uniqueLocIdxSlice, 1 ) );

numUniqueIdxs = transpose( sum( uniqueLocIdxSlice, 1 ) );
invNumUniqueIdxs = 1 ./ numUniqueIdxs;

Mxavg_inslice_xy = invNumUniqueIdxs .* Mxsum_inslice_xy;
Myavg_inslice_xy = invNumUniqueIdxs .* Mysum_inslice_xy;

Mxytargph_inslice = atan2( Myavg_inslice_xy, Mxavg_inslice_xy );
targFAphase_inslice = exp( 1j * Mxytargph_inslice );

thetatarg_inslice = targFA_inslice .* targFAphase_inslice( uniqueIdxToPos );

compFAthru = zeros( numPosSubj, 1 );
compFAtarg = zeros( numPosSubj, 1 );

compFAthru( subjNSlicesLocIdx ) = theta_outslice;
compFAtarg( subjNSlicesLocIdx ) = thetatarg_outslice;

compFAthru( subjSlicesLocIdx ) = theta_inslice;
compFAtarg( subjSlicesLocIdx ) = thetatarg_inslice;

magNRMSE_compFA = norm( compFAthru - compFAtarg, 2 ) / norm( compFAtarg, 2 );

% calculate magNRMSE_FA
magNRMSE_FA = norm( FAthru - FAtarg, 2 ) / norm( FAtarg, 2 );

for ss = 1:numSlices
    
    % find slice bounds
    inSliceIdx = find(...
        ( val.pos( :, 3 ) >= ( sliceBounds( ss, 1 ) - di_tol ) ) &...
        ( val.pos( :, 3 ) <= ( sliceBounds( ss, 2 ) + di_tol ) ) );

    [ ~, sliceCenterPosIdx ] = ...
        min( vecnorm( val.pos( inSliceIdx, : ) - [ 0, 0, sliceLocation( ss ) ], 2, 2 ) );

    [ ~, sliceCenterIdx ] = min( abs( val.z - sliceLocation( ss ) ) );

    targValues = Mtarg( inSliceIdx( sliceCenterPosIdx ), : );
    targFA = real( acosd( targValues( 1, 3 ) ) );

    % increase slice width to compute FWHM and transitionWidth
    sliceThickness = abs( diff( sliceBounds( ss, : ) ) );
    sliceBoundsExtend = [...
        sliceBounds( ss, 1 ) - sliceBoundExtendAmount * sliceThickness,...
        sliceBounds( ss, 2 ) + sliceBoundExtendAmount * sliceThickness,...
        ];

    sliceExtendIdx = find(...
        ( val.pos( :, 3 ) >= ( sliceBoundsExtend( 1 ) - di_tol ) ) &...
        ( val.pos( :, 3 ) <= ( sliceBoundsExtend( 2 ) + di_tol ) ) );

    sliceExtend = [...
        find( ( val.z >= ( sliceBoundsExtend( 1 ) - di_tol ) ), 1, 'first' ),...
        find( ( val.z <= ( sliceBoundsExtend( 2 ) + di_tol ) ), 1, 'last' )...
        ];

    % sliceExtendFull = sliceExtend(1) : 1 : sliceExtend(end);
    sliceExtendBounds = [ val.z( sliceExtend( 1 ) ), val.z( sliceExtend( end ) ) ];

    sliceExtendLeft = sliceExtend( 1 ):sliceCenterIdx;
    sliceExtendLeftBounds = [ val.z( sliceExtendLeft( 1 ) ), val.z( sliceExtendLeft( end ) ) ];

    sliceExtendRight = sliceCenterIdx:sliceExtend( 2 );
    sliceExtendRightBounds = [ val.z( sliceExtendRight( 1 ) ), val.z( sliceExtendRight( end ) ) ];
    
    % determine unique XY positions to use for metrics
    MxThruSliceExtend = repmat( MThru( sliceExtendIdx, 1 ), [ 1, val.numUniqueXYPos ] );
    MyThruSliceExtend = repmat( MThru( sliceExtendIdx, 2 ), [ 1, val.numUniqueXYPos ] );
    FAThruSliceExtend = repmat( FAthru( sliceExtendIdx ), [ 1, val.numUniqueXYPos ] );
    FAtargSliceExtend = repmat( FAtarg( sliceExtendIdx ), [ 1, val.numUniqueXYPos ] );
    zThruSliceExtend = repmat( val.pos( sliceExtendIdx, 3 ), [ 1, val.numUniqueXYPos ] );

    numInSlice = size( sliceExtendIdx, 1 );
    uniqueIdxToPos = val.uniqueXYPosToXYPos( sliceExtendIdx, : );
    uniquePosXYIdxsMat = logical( sparse(...
        uint32(transpose(1:numInSlice)),...
        uniqueIdxToPos,...
        ones( numInSlice, 1 ),...
        numInSlice, val.numUniqueXYPos ) );

    uniqueIdxCts = full( sum( uniquePosXYIdxsMat, 1 ) );
    maxUniqueIdxCts = max( uniqueIdxCts );
    uniqueIdxsToUse = uniqueIdxCts == maxUniqueIdxCts;
    numUniqueIdxsToUse = sum( uniqueIdxsToUse );
    
    uniquePosXYIdxsMatUpdate = uniquePosXYIdxsMat;
    uniquePosXYIdxsMatUpdate( :, ~uniqueIdxsToUse ) = false;
    
    xySamplePtMxVals{ ss } = reshape( MxThruSliceExtend( uniquePosXYIdxsMatUpdate ),...
       [ maxUniqueIdxCts, numUniqueIdxsToUse ] );
    xySamplePtMyVals{ ss } = reshape( MyThruSliceExtend( uniquePosXYIdxsMatUpdate ),...
       [ maxUniqueIdxCts, numUniqueIdxsToUse ] );
    xySamplePtFAVals{ ss } = reshape( FAThruSliceExtend( uniquePosXYIdxsMatUpdate ),...
       [ maxUniqueIdxCts, numUniqueIdxsToUse ] );
    xySamplePtFAtargVals{ ss } = reshape( FAtargSliceExtend( uniquePosXYIdxsMatUpdate ),...
       [ maxUniqueIdxCts, numUniqueIdxsToUse ] );
    xySamplePtzVals{ ss } = reshape( zThruSliceExtend( uniquePosXYIdxsMatUpdate ),...
       [ maxUniqueIdxCts, numUniqueIdxsToUse ] );

    [ xySamplePtzVals{ ss }, xySamplePtzValsSortIdx ] = sort( xySamplePtzVals{ ss }, 1 );
    xySamplePtzValsSortLinIdx = sub2ind( [ maxUniqueIdxCts, numUniqueIdxsToUse ],...
        xySamplePtzValsSortIdx, repmat( 1:numUniqueIdxsToUse, [ maxUniqueIdxCts, 1 ] ) );
    xySamplePtFAVals{ ss }( : ) = xySamplePtFAVals{ ss }( xySamplePtzValsSortLinIdx );
    xySamplePtMxVals{ ss }( : ) = xySamplePtMxVals{ ss }( xySamplePtzValsSortLinIdx );
    xySamplePtMyVals{ ss }( : ) = xySamplePtMyVals{ ss }( xySamplePtzValsSortLinIdx );

    xySamplePtPhVals{ ss } = unwrap( atan2( xySamplePtMyVals{ ss }, xySamplePtMxVals{ ss } ), [], 1 );
    
    transitionWidthUnique = zeros( numUniqueIdxsToUse, 1 );
    FWHMUnique = zeros( numUniqueIdxsToUse, 1 );
    phaseChangeAcrossSliceUnique = zeros( numUniqueIdxsToUse, 1 );
    inSliceRippleUnique = zeros( numUniqueIdxsToUse, 1 );
    outSliceRippleUnique = zeros( numUniqueIdxsToUse, 1 );
    inSliceRipplePercentUnique = zeros( numUniqueIdxsToUse, 1 );
    outSliceRipplePercentUnique = zeros( numUniqueIdxsToUse, 1 );

    for uu = 1:numUniqueIdxsToUse

        % create phase sampling functions with unwrapping
        phThruSampleFun = makeWaveformSampleFunction( xySamplePtzVals{ ss }( :, 1 ), xySamplePtPhVals{ ss }( :, uu ), interpMethod, extrap );

        % Create FA targ sampling functions
        FAtargSliceSampleFunFull = makeWaveformSampleFunction(...
            xySamplePtzVals{ ss }( :, uu ), xySamplePtFAVals{ ss }( :, uu ),...
            interpMethod, extrap );

        zValsIdxLeft = xySamplePtzVals{ ss }( :, uu ) <= ( val.sliceLocation + di_tol );
        zValsIdxRight = xySamplePtzVals{ ss }( :, uu ) >= ( val.sliceLocation - di_tol );

        maxFAUniqueLeft = min( max( xySamplePtFAVals{ ss }( zValsIdxLeft, uu ) ), targFA );
        maxFAUniqueRight = min( max( xySamplePtFAVals{ ss }( zValsIdxRight, uu ) ), targFA );
        minFAUniqueLeft = max( min( xySamplePtFAVals{ ss }( zValsIdxLeft, uu ) ) + FApercentage * targFA, TBWLow * targFA );
        minFAUniqueRight = max( min( xySamplePtFAVals{ ss }( zValsIdxRight, uu ) ) + FApercentage * targFA, TBWLow * targFA );

        halfValLeft = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - FWHMVal * targFA, sliceExtendLeftBounds, sliceExtendBounds, sliceBounds, 'left' );
        highValLeft = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - TBWHigh * maxFAUniqueLeft, sliceExtendLeftBounds, sliceExtendBounds, sliceBounds, 'left' );
        lowValLeft = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - minFAUniqueLeft, sliceExtendLeftBounds, sliceExtendBounds, sliceBounds, 'left' );

        halfValRight = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - FWHMVal * targFA, sliceExtendRightBounds, sliceExtendBounds, sliceBounds, 'right' );
        highValRight = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - TBWHigh * maxFAUniqueRight, sliceExtendRightBounds, sliceExtendBounds, sliceBounds, 'right' );
        lowValRight = findRootsInterpolation( @( z ) FAtargSliceSampleFunFull( z ) - minFAUniqueRight, sliceExtendRightBounds, sliceExtendBounds, sliceBounds, 'right' );

        transitionWidthUnique( uu ) = mean( [ abs( highValLeft - lowValLeft ); abs( highValRight - lowValRight ) ] );
        FWHMUnique( uu ) = mean( halfValRight - halfValLeft );

        phaseChangeAcrossSliceUnique( uu ) = mean( abs( diag( phThruSampleFun( halfValLeft ) ) - diag( phThruSampleFun( halfValRight ) ) ) );

        % find points that define slices
        sliceInIdx = sort( [...
            find( ( xySamplePtzVals{ ss }( :, uu ) >= ( highValLeft - di_tol ) ), 1, 'first' ),...
            find( ( xySamplePtzVals{ ss }( :, uu ) <= ( highValRight + di_tol ) ), 1, 'last' )...
            ] );
        sliceOutIdx = sort( [...
            find( ( xySamplePtzVals{ ss }( :, uu ) >= ( lowValLeft - di_tol ) ), 1, 'first' ) - 1,...
            find( ( xySamplePtzVals{ ss }( :, uu ) <= ( lowValRight + di_tol ) ), 1, 'last' ) + 1 ...
            ] );

        FAThruSlice = xySamplePtFAVals{ ss }( sliceInIdx(1):sliceInIdx(2), uu );
        FAThruOutOfSlice = [...
            xySamplePtFAVals{ ss }( 1:( sliceOutIdx( 1 ) ), uu );...
            xySamplePtFAVals{ ss }( ( sliceOutIdx( 2 ) ):end, uu ) ];

        inSliceRippleUnique( uu ) = max( FAThruSlice ) - min( FAThruSlice );
        outSliceRippleUnique( uu ) = max( FAThruOutOfSlice ) - min( FAThruOutOfSlice );

        inSliceRipplePercentUnique( uu ) = inSliceRippleUnique( uu ) / targFA;
        outSliceRipplePercentUnique( uu ) = outSliceRippleUnique( uu ) / targFA;
    end

    transitionWidth( ss ) = mean( transitionWidthUnique );
    FWHM( ss ) = mean( FWHMUnique );
    phaseChangeAcrossSlice( ss ) = mean( phaseChangeAcrossSliceUnique );

    inSliceRipple( ss ) = mean( inSliceRippleUnique );
    outSliceRipple( ss ) = mean( outSliceRippleUnique );

    inSliceRipplePercent( ss ) = mean( inSliceRipplePercentUnique );
    outSliceRipplePercent( ss ) = mean( outSliceRipplePercentUnique );

end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function r = findRootsInterpolation( fun, boundsbias, boundstotal, sliceBounds, shift )

% sample function with a certain amount of points and find the places where
% value is near 0
zeroTol = 1e-2;
numSampPtsBias = 1e6;
samplePtsBoundsBias = linspace( boundsbias( 1 ), boundsbias( 2 ), numSampPtsBias );

% try to determine the points closest to 0

zeroThreshBias = abs( fun( samplePtsBoundsBias ) ) < zeroTol;
% [ ~, zeroApproxBiasIdxs ] = max( zeroThresh, [], 1 );

if matches( shift, 'left' )
    sideIdx = 1;
elseif matches( shift, 'right' )
    sideIdx = 2;
else
    error( "Unknown shift." )
end

if any( zeroThreshBias, 1 )
    
    threshIdxs = find( zeroThreshBias );
    [ ~, threshIdx ] = min( abs( samplePtsBoundsBias( threshIdxs ) - sliceBounds( sideIdx ) ) );
    r = samplePtsBoundsBias( threshIdxs( threshIdx ) );

else % there isn't a root in that domain.... make the domain bigger

    numSampPtsTotal = 10 * numSampPtsBias;
    samplePtsBoundsTotal = linspace( boundstotal( 1 ), boundstotal( 2 ), numSampPtsTotal );

    % find points closest to 0
    
    absFunSampleTotal = abs( fun( samplePtsBoundsTotal ) );
    zeroThreshTotal = absFunSampleTotal < zeroTol;
    % [ zeroApproxTotalIdxs ] = find( absFunSampleTotal < zeroTol );

    if any( zeroThreshTotal, 1 )
        threshIdxs = find( zeroThreshTotal );
        [ ~, threshIdx ] = min( abs( samplePtsBoundsTotal( threshIdxs ) - sliceBounds( sideIdx ) ) );
        r = samplePtsBoundsTotal( threshIdxs( threshIdx ) );
    else
        [ ~, minApproxIdx ] = min( absFunSampleTotal, [], 1 );
        r = samplePtsBoundsTotal( minApproxIdx );
    end

end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function r = findRootsChebApprox( fun, boundsbias, boundstotal, shift ) %#ok
tol = 1e-12;
maxCheb = 1000;

% Project function into chebyshev basis
c = chebCoeffFromFun( fun, boundsbias, tol );

% Find roots of the associated Chebyshev polynomial
rs = chebRoots( c, boundsbias );

% select the root that corresponds to the direction
if ~isempty( rs )
    if matches( shift, 'left' )
        r = rs( 1 );
    elseif matches( shift, 'right' )
        r = rs( end );
    else
        error( "Unknown shift." )
    end
else
    cfull = chebCoeffFromFun( fun, boundstotal, tol, maxCheb );
    rs = chebRoots( cfull, boundstotal );
    if ~isempty( rs )
        if matches( shift, 'right' )
            r = rs( 1 );
        elseif matches( shift, 'left' )
            r = rs( end );
        else
            error( "Unknown shift." )
        end
    else
        [ ~, rs ] = chebMinMax( cfull, boundstotal );
        r = rs( end );
    end
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function x0 = findRootsfzero( fun, bounds, shift, fzeroOption ) %# ok

maxIter = 10;

if matches( shift, 'left' )
    shiftIdx = [ true, false ];
elseif matches( shift, 'right' )
    shiftIdx = [ false, true ];
else
    error( "Unknown shift." )
end

diffTol = 1e-8;
weightAvg = 0.9;

x0p = fzero(...
    fun, bounds, fzeroOption );
newBoundsFn = @( x0p ) sort( [ bounds( shiftIdx ), ( x0p * weightAvg + ( 1 - weightAvg ) * bounds( ~shiftIdx ) ) ] );

ii = 0;
while true
    ii = ii + 1;
    newBounds = newBoundsFn( x0p );
    if sign( fun( newBounds( 1 ) ) * fun( newBounds( 2 ) ) )  ~= -1
        ishappy = false;
        break;
    end

    x0 = fzero(...
        fun, newBounds, fzeroOption );

    if abs( x0 - x0p ) < diffTol
        ishappy = true;
        break;
    else
        x0p = x0;
    end

    if ii >= maxIter
        ishappy = false;
        break;
    end

end

if ~ishappy
    x0 = fzero(...
        fun, bounds, fzeroOption );
    warning( "Couldn't converge to fzero point." )
end
end
% ----------------------------------------------------------------------- %