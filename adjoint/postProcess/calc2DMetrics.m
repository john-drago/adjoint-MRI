function [ val ] = calc2DMetrics( val, pulse, fields )

si = val.si;

%% Assign workers
% pp = gcp('nocreate');
% numcores = feature('numcores');
% numWorkers = 8;
% 
% if isempty(pp)
%     parpool( min( [ numWorkers, numcores ] ) );
% else
%     if pp.NumWorkers ~= numWorkers
%         delete( pp );
%         parpool( min( [ numWorkers, numcores ] ) );
%     end
% end

%% Need to simulate through slice behavior
thruplane_ss_di = 0.1e-3;
thruplane_dxdy_di = 20e-3;

inplane_ss_di = 1e-3;
inplane_dxdy_di = 1e-3;

val.MThru = cell( 1, val.numSubj );
val.MtargThru = cell( 1, val.numSubj );
val.FAThru = cell( 1, val.numSubj );
val.FAtargThru = cell( 1, val.numSubj );
val.zThru = cell( 1, val.numSubj );

val.calcThruMetrics = true( pulse.numSlices, val.numSubj );

val.inSliceRipple = zeros( pulse.numSlices, val.numSubj );
val.outOfSliceRipple = zeros( pulse.numSlices, val.numSubj );
val.inSliceRipplePercent = zeros( pulse.numSlices, val.numSubj );
val.outOfSliceRipplePercent = zeros( pulse.numSlices, val.numSubj );
val.transitionWidth = zeros( pulse.numSlices, val.numSubj );
val.phaseChangeAcrossSlice = zeros( pulse.numSlices, val.numSubj );
val.FWHM = zeros( pulse.numSlices, val.numSubj );
val.magNRMSEThruPlane_FA = zeros( 1, val.numSubj );
val.magNRMSEThruPlane_compFA = zeros( 1, val.numSubj );
val.xySamplePtzVals = cell( pulse.numSlices, val.numSubj );
val.xySamplePtFAtargVals = cell( pulse.numSlices, val.numSubj );
val.xySamplePtFAVals = cell( pulse.numSlices, val.numSubj );
val.xySamplePtPhVals = cell( pulse.numSlices, val.numSubj );
val.xySamplePtMxVals = cell( pulse.numSlices, val.numSubj );
val.xySamplePtMyVals = cell( pulse.numSlices, val.numSubj );

%% Pre-allocate

fieldsX = fields.X;
fieldsY = fields.Y;
fieldsZ = fields.Z;
fieldsb1p = fields.b1p;
fieldsbz = fields.bz;
pulseshimarray = pulse.shimarray;
fieldsM0 = fields.M0;
fieldsMtarg = fields.Mtarg;

%% Interpolate into the points of interest
interpMtarg = 'nearest';
di_tol = thruplane_ss_di / 10;
maskInterpThresh = 0.50 - eps;

for ss = 1:val.numSubj

    % fprintf( "\nCalc2DMetrics Thru Slice %s\n", val.subjIden( ss )  );

    ssi = si( ss );

    % Determine the z points within the roi
    zOpt = unique( fields.Z( fields.opt_roi( :, :, :, ssi ) ) );
    zThru = unique( [ 0, zOpt( 1 ) : thruplane_ss_di : zOpt( end ) ] );
    xUnique = unique( fields.X( fields.opt_roi( :, :, :, ssi ) ) );
    yUnique = unique( fields.Y( fields.opt_roi( :, :, :, ssi ) ) );
    xThru = unique( [ 0, xUnique( 1 ) : thruplane_dxdy_di : xUnique( end ) ] );
    yThru = unique( [ 0, yUnique( 1 ) : thruplane_dxdy_di : yUnique( end ) ] );
    % xThru = 0;
    % yThru = 0;
    val.zThru{ ss } = zThru;

    [ XThru, YThru, ZThru ] = ndgrid( xThru, yThru, zThru );

    % Interpolate optROI
    optROIInterpolate = Interp3D(...
        fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, :, ssi ) ) ),...
        XThru, YThru, ZThru ) > maskInterpThresh;
    numThruROI = sum( optROIInterpolate, "all" );

    % Interpolate valROI
    % valROIInterpolate = Interp3D(...
    %     fields.X, fields.Y, fields.Z, double( squeeze( fields.val_roi( :, :, :, ssi ) ) ),...
    %     val.X, val.Y, val.Z ) > maskInterpThresh;

    % Interpolate b1p
    b1pThru = zeros( numThruROI, val.numXYCoils );
    % parfor nn = 1:val.numXYCoils
    for nn = 1:val.numXYCoils

        % b1pROI_nn = squeeze(...
        %     Interp3D(...
        %     fields.X, fields.Y, fields.Z, squeeze( fields.b1p( nn, :, :, :, ssi ) ),...
        %     XThru, YThru, ZThru ) );

        b1pROI_nn = ...
            Interp3D(...
            fieldsX, fieldsY, fieldsZ, squeeze( fieldsb1p( nn, :, :, :, ssi ) ),...
            XThru, YThru, ZThru );
        b1pThru( :, nn ) = squeeze( b1pROI_nn( optROIInterpolate ) );
    end

    % Interpolate db0
    
    % db0ROI = squeeze(...
    %     Interp3D(...
    %     fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, :, si(ss) ) ),...
    %         XThru, YThru, ZThru ) );

    db0ROI = ...
        Interp3D(...
        fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, :, si(ss) ) ),...
            XThru, YThru, ZThru );

    db0Thru = squeeze( db0ROI( optROIInterpolate ) );

    % Interpolate bz
    if val.numZCoils > 0

        bzThru = zeros( numThruROI, val.numZCoils );

        for nn = 1:val.numZCoils
        % parfor nn = 1:val.numZCoils
            
            bzROI = []; %#ok
            
            if pulseshimarray
                % bzROI = squeeze(...
                %     Interp3D(...
                %     fields.X, fields.Y, fields.Z, squeeze( fields.bz( nn, :, :, : ) ),...
                %     XThru, YThru, ZThru ) );

                bzROI = ...
                    Interp3D(...
                    fieldsX, fieldsY, fieldsZ, squeeze( fieldsbz( nn, :, :, : ) ),...
                    XThru, YThru, ZThru );

            else
                error( "Need to assign to pulse.shimarray" );
            end

            bzThru( :, nn ) = squeeze( bzROI( optROIInterpolate ) );
        end
    else
        bzThru = zeros( numThruROI, 1 );
    end

    % Interpolate Mtarg
    MtargThru = zeros( numThruROI, 3 );

    for dd = 1:3

        % MtargThruInitInterp_dd = squeeze( Interp3D(...
        %     fields.X, fields.Y, fields.Z, squeeze( fields.Mtarg( :, :, :, dd ) ),...
        %     XThru, YThru, ZThru, interpMtarg ) );

        MtargThruInitInterp_dd = Interp3D(...
            fields.X, fields.Y, fields.Z, squeeze( fields.Mtarg( :, :, :, dd ) ),...
            XThru, YThru, ZThru, interpMtarg );

        MtargThru_dd = zeros( size( MtargThruInitInterp_dd ) );

        for ll = 1 : val.numSlices

            if ( ll == 1 ) && ( ll == val.numSlices )
                lowBoundIdx = 1;
                highBoundIdx = length( zThru );
            elseif ll == 1
                lowBoundIdx = 1;
                highBoundIdx = find( zThru <= ( val.sliceBounds( ll+1, 1 ) - di_tol ), 1, 'last' );
            elseif ll == val.numSlices
                lowBoundIdx = find( zThru >= ( val.sliceBounds( ll-1, 2 ) + di_tol ), 1, 'first' );
                highBoundIdx = length( zThru );
            end

            sliceIdxStart = find( zThru >= ( val.sliceBounds( ll, 1 ) - di_tol ), 1, 'first' );
            sliceIdxEnd = find( zThru <= ( val.sliceBounds( ll, 2 ) + di_tol ), 1, 'last' );

            if ( sliceIdxEnd == length( zThru ) ) && isempty( sliceIdxStart )
                sliceIdxStart = highBoundIdx;
                sliceIdxEnd = highBoundIdx;
            elseif ( sliceIdxStart == 1 ) && isempty( sliceIdxEnd )
                sliceIdxStart = lowBoundIdx;
                sliceIdxEnd = lowBoundIdx;
            end

            if sliceIdxStart == sliceIdxEnd
                val.calcThruMetrics( ll, ss ) = false;
            end

            sliceMiddleIdx = ceil( mean( unique( [ sliceIdxStart, sliceIdxEnd ] ) ) );
            if sliceIdxStart ~= lowBoundIdx
                sliceLeftIdx = floor( mean( unique( [ lowBoundIdx, sliceIdxStart-1 ] ) ) );
            else
                sliceLeftIdx = lowBoundIdx;
            end
            if sliceIdxEnd ~= highBoundIdx
                sliceRightIdx = ceil( mean( [ sliceIdxEnd+1, highBoundIdx ] ) );
            else
                sliceRightIdx = highBoundIdx;
            end

            outsideSliceValLeft = MtargThruInitInterp_dd( :, :, sliceLeftIdx );
            outsideSliceValRight = MtargThruInitInterp_dd( :, :, sliceRightIdx );
            sliceVal = MtargThruInitInterp_dd( :, :, sliceMiddleIdx );

            if ~isempty( outsideSliceValLeft )
                MtargThru_dd( :, :, lowBoundIdx:(sliceIdxStart-1) ) = repmat( outsideSliceValLeft,...
                    [ 1, 1, sliceIdxStart-lowBoundIdx ] );
            end
            if ~isempty( sliceVal )
                MtargThru_dd( :, :, sliceIdxStart:sliceIdxEnd ) = repmat( sliceVal,...
                    [ 1, 1, sliceIdxEnd-sliceIdxStart+1 ] );
            end
            if ~isempty( outsideSliceValRight )
                MtargThru_dd( :, :, (sliceIdxEnd+1):highBoundIdx ) = repmat( outsideSliceValRight,...
                    [ 1, 1, highBoundIdx-sliceIdxEnd ] );
            end

        end

        MtargThru( :, dd ) = squeeze( MtargThru_dd( optROIInterpolate ) );
    end

    FAtargThru = real( acosd( MtargThru( :, 3 ) ) );

    % Interpolate M0
    M0Thru = zeros( numThruROI, 3 );
    % parfor dd = 1:3
    for dd = 1:3
        M0_dd = squeeze( Interp3D(...
            fieldsX, fieldsY, fieldsZ, squeeze( fieldsM0( :, :, :, dd ) ),...
            XThru, YThru, ZThru, interpMtarg ) );
        M0Thru( :, dd ) = squeeze( M0_dd( optROIInterpolate ) );
    end

    % get the pos
    PosThru = [...
        squeeze( XThru( optROIInterpolate ) ),...
        squeeze( YThru( optROIInterpolate ) ),...
        squeeze( ZThru( optROIInterpolate ) ) ];

    % Get the waveform
    wv = val.generateWaveforms( val.pOpt, val );

    % Create new opt struct
    thruOpt = struct;
    thruOpt.b1preal = real( b1pThru );
    thruOpt.b1pimag = imag( b1pThru );
    thruOpt.bzsens = bzThru;
    thruOpt.db0 = db0Thru;
    thruOpt.pos = PosThru;
    thruOpt.M0 = M0Thru;
    thruOpt.Mtarg = MtargThru;
    thruOpt.gyro = val.gyro;

    thruOpt.x = xThru;
    thruOpt.y = yThru;
    thruOpt.z = zThru;

    thruOpt.numZCoils = val.numZCoils;
    thruOpt.numXYCoils = val.numXYCoils;
    thruOpt.numSubj = val.numSubj;
    thruOpt.convertMBackToLarmor = false;

    thruOpt.dwxyvec = val.dwxyvec;
    thruOpt.numTimePoints = val.numTimePoints;
    thruOpt.numPos = numThruROI;

    thruOpt.sliceLocation = val.sliceLocation;
    thruOpt.sliceBounds = val.sliceBounds;
    thruOpt.numSlices = val.numSlices;

    [ thruOpt, pulse ] = processSliceSelectivePulse( thruOpt, pulse );
    thruOpt.numSubj = 1;
    thruOpt.cum_numPosSubj = size( thruOpt.pos, 1 );
    thruOpt = getSliceIdxs( thruOpt );
    [ thruOpt.uniqueXYPos, thruOpt.XYPosToUniqueXYPos, thruOpt.uniqueXYPosToXYPos, thruOpt.numUniqueXYPos  ] = ...
        getUniqueXYPos( thruOpt.pos );

    saveIntermediate = false;
    val.MThru{ ss } = runAdjointForwardModelCPU( wv, thruOpt, saveIntermediate );
    val.FAThru{ ss } = real( acosd( val.MThru{ ss }( :, 3 ) ) );
    val.MtargThru{ ss } = MtargThru;
    val.FAtargThru{ ss } = FAtargThru;

    if all( val.calcThruMetrics( :, ss ) )
        % calculate metrics
        [ val.inSliceRipple( :, ss ), val.outOfSliceRipple( :, ss ), val.inSliceRipplePercent( :, ss ), val.outOfSliceRipplePercent( :, ss ),...
            val.transitionWidth( :, ss ), val.FWHM( :, ss ), val.phaseChangeAcrossSlice( :, ss ), val.magNRMSEThruPlane_FA( ss ), val.magNRMSEThruPlane_compFA( ss ),...
            val.xySamplePtMxVals( :, ss ), val.xySamplePtMyVals( :, ss ),...
            val.xySamplePtFAVals( :, ss ), val.xySamplePtPhVals( :, ss ),...
            val.xySamplePtzVals( :, ss ), val.xySamplePtFAtargVals( :, ss ) ] =...
            calculateThruSliceMetrics(...
            val.MThru{ ss }, val.MtargThru{ ss }, val.FAThru{ ss }, val.FAtargThru{ ss }, thruOpt );
    end

end

%% Need to simulate in slice behavior
maskInterpThresh = 0.50 - eps;

val.MvecInPlane = cell( pulse.numSlices, val.numSubj );
val.FAvecInPlane = cell( pulse.numSlices, val.numSubj );
val.MInPlane = cell( pulse.numSlices, val.numSubj );
val.FAInPlane = cell( pulse.numSlices, val.numSubj );
val.optroiInPlane = cell( pulse.numSlices, val.numSubj );
val.zInPlane = cell( pulse.numSlices, val.numSubj );

val.calcInPlaneMetrics = true( pulse.numSlices, val.numSubj );

val.magNRMSEInPlane_Mxy = zeros( pulse.numSlices, val.numSubj );
val.NRMSEInPlane_Mxy = zeros( pulse.numSlices, val.numSubj );
val.magNRMSEInPlane_FA = zeros( pulse.numSlices, val.numSubj );
val.meanFAInPlane = zeros( pulse.numSlices, val.numSubj );
val.stdFAInPlane = zeros( pulse.numSlices, val.numSubj );
val.nMaxMinInPlane = zeros( pulse.numSlices, val.numSubj );
val.magCoeffVarInPlane = zeros( pulse.numSlices, val.numSubj );

val.magNRMSEInPlaneCentralSlice_Mxy = zeros( pulse.numSlices, val.numSubj );
val.NRMSEInPlaneCentralSlice_Mxy = zeros( pulse.numSlices, val.numSubj );
val.magNRMSEInPlaneCentralSlice_FA = zeros( pulse.numSlices, val.numSubj );
val.meanFAInPlaneCentralSlice = zeros( pulse.numSlices, val.numSubj );
val.stdFAInPlaneCentralSlice = zeros( pulse.numSlices, val.numSubj );
val.nMaxMinInPlaneCentralSlice = zeros( pulse.numSlices, val.numSubj );
val.magCoeffVarInPlaneCentralSlice = zeros( pulse.numSlices, val.numSubj );

numZCoordInSlice = floor( pulse.sliceThickness / inplane_ss_di );

for ll = 1:pulse.numSlices

    % Determine the z location in which to perform roi interpolation
    sliceBounds = pulse.sliceBounds( ll, : );
    sliceDiscretization = ( diff( sliceBounds ) - ( numZCoordInSlice - 1 ) * inplane_ss_di ) / 2 + ( inplane_ss_di * ( 0 : ( numZCoordInSlice - 1 ) ) );
    zInPlane = unique( [ sliceBounds, sliceBounds( ll, 1 ) + sliceDiscretization, pulse.sliceLocation( ll ) ] );

    xInPlane = val.x( 1 ) : inplane_dxdy_di : val.x( end );
    yInPlane = val.y( 1 ) : inplane_dxdy_di : val.y( end );

    [ XInPlane, YInPlane, ZInPlane ] = ndgrid( xInPlane, yInPlane, zInPlane );

    for ss = 1:val.numSubj

        % fprintf( "\nCalc2DMetrics In Slice %s\n", val.subjIden( ss )  );

        optROIInPlane = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, :, si(ss) ) ) ),...
            XInPlane, YInPlane, ZInPlane ) > maskInterpThresh;

        val.optroiInPlane{ ll, ss } = optROIInPlane;

        idxVec = uint32( find( optROIInPlane( : ) ) );
        [ ri, ci, di ] = ind2sub( size( optROIInPlane ), idxVec );

        % Interpolate b1p
        b1pInPlane = zeros( length( idxVec ), val.numXYCoils );
        % parfor nn = 1:val.numXYCoils
        for nn = 1:val.numXYCoils
            b1pArr = squeeze(...
                Interp3D(...
                fieldsX, fieldsY, fieldsZ, squeeze( fieldsb1p( nn, :, :, :, ssi ) ),...
                XInPlane, YInPlane, ZInPlane ) );
            b1pInPlane( :, nn ) = b1pArr( optROIInPlane );
        end

        % Interpolate db0
        db0Arr = squeeze(...
            Interp3D(...
            fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, :, si(ss) ) ),...
            XInPlane, YInPlane, ZInPlane ) );

        db0InPlane = db0Arr( optROIInPlane );

        % Interpolate bz
        if val.numZCoils > 0
            bzInPlane = zeros( length( idxVec ), val.numZCoils );

            % parfor nn = 1:val.numZCoils
            for nn = 1:val.numZCoils

                bzArr = []; %#ok

                if pulseshimarray
                    bzArr = squeeze(...
                        Interp3D(...
                        fieldsX, fieldsY, fieldsZ, squeeze( fieldsbz( nn, :, :, : ) ),...
                        XInPlane, YInPlane, ZInPlane ) );
                else
                    error( "Need to assign to pulse.shimarray" );
                end

                bzInPlane( :, nn ) = bzArr( optROIInPlane );
            end
        else
            bzInPlane = zeros( length( idxVec ), 1 );
        end

        % Interpolate Mtarg
        MtargInPlane = zeros( length( idxVec ), 3 );
        % parfor dd = 1:3
        for dd = 1:3
            MtargArr = squeeze(...
                Interp3D(...
                fieldsX, fieldsY, fieldsZ, squeeze( fieldsMtarg( :, :, :, dd ) ),...
                XInPlane, YInPlane, ZInPlane ) );
            MtargInPlane( :, dd ) = MtargArr( optROIInPlane );
        end
        FAtargInPlane = real( acosd( MtargInPlane( :, 3 ) ) );

        % Interpolate M0
        M0InPlane = zeros( length( idxVec ), 3 );
        % parfor dd = 1:3
        for dd = 1:3
            M0Arr = squeeze(...
                Interp3D(...
                fieldsX, fieldsY, fieldsZ, squeeze( fieldsM0( :, :, :, dd ) ),...
                XInPlane, YInPlane, ZInPlane ) );

            M0InPlane( :, dd ) = M0Arr( optROIInPlane );
        end

        % get the pos
        PosInPlane = [ XInPlane( optROIInPlane ), YInPlane( optROIInPlane ), ZInPlane( optROIInPlane ) ];

        % Get the waveform
        wv = val.generateWaveforms( val.pOpt, val );

        % Create new opt struct
        inplaneOpt = struct;
        inplaneOpt.b1preal = real( b1pInPlane );
        inplaneOpt.b1pimag = imag( b1pInPlane );
        inplaneOpt.bzsens = bzInPlane;
        inplaneOpt.db0 = db0InPlane;
        inplaneOpt.pos = PosInPlane;
        inplaneOpt.M0 = M0InPlane;
        inplaneOpt.Mtarg = MtargInPlane;
        inplaneOpt.gyro = val.gyro;

        inplaneOpt.numZCoils = val.numZCoils;
        inplaneOpt.numXYCoils = val.numXYCoils;
        inplaneOpt.numSubj = val.numSubj;
        
        inplaneOpt.dwxyvec = val.dwxyvec;
        inplaneOpt.numTimePoints = val.numTimePoints;
        inplaneOpt.numPos = length( zThru );

        saveIntermediate = false;
        val.MvecInPlane{ ll, ss } = runAdjointForwardModelCPU( wv, inplaneOpt, saveIntermediate );
        val.FAvecInPlane{ ll, ss } = real( acosd( val.MvecInPlane{ ll, ss }( :, 3 ) ) );
        
        val.MInPlane{ ll, ss } = zeros( [ size( optROIInPlane ), 3 ] );
        val.FAInPlane{ ll, ss } = zeros( [ size( optROIInPlane ) ] );
        
        val.FAInPlane{ ll, ss }( idxVec ) = val.FAvecInPlane{ ll, ss };

        for dd = 1:3
            li = sub2ind( [ size( optROIInPlane ), 3 ], ri, ci, di, dd*ones(size(ri)) );
            val.MInPlane{ ll, ss }( li ) = val.MvecInPlane{ ll, ss }( :, dd );
        end

        val.zInPlane{ ll, ss } = zInPlane;

        % calculate metrics for entire inplane slice
        if ~all( ~optROIInPlane, 'all' )
            [ val.magNRMSEInPlane_Mxy( ll, ss ), val.NRMSEInPlane_Mxy( ll, ss ), val.magNRMSEInPlane_FA( ll, ss ),...
                val.meanFAInPlane( ll, ss ), val.stdFAInPlane( ll, ss ), val.nMaxMinInPlane( ll, ss ), val.magCoeffVarInPlane( ll, ss ) ] = ...
                calculateFAMtargMetrics( val.MvecInPlane{ ll, ss }, val.FAvecInPlane{ ll, ss }, MtargInPlane, FAtargInPlane );
        else
            val.calcInPlaneMetrics( ll, ss ) = false;
        end

        sliceLocIdxs = abs( PosInPlane( :, 3 ) - pulse.sliceLocation ) < di_tol;

        MvecInPlaneCentralSlice = val.MvecInPlane{ ll, ss }( sliceLocIdxs, : );
        FAvecInPlaneCentralSlice = val.FAvecInPlane{ ll, ss }( sliceLocIdxs, : );
        MtargInPlaneCentralSlice = MtargInPlane( sliceLocIdxs, : );
        FAtargInPlaneCentralSlice = FAtargInPlane( sliceLocIdxs, : );

        % calculate metrics for central inplane slice
        if ~all( ~sliceLocIdxs, 'all' )
            [ val.magNRMSEInPlaneCentralSlice_Mxy( ll, ss ), val.NRMSEInPlaneCentralSlice_Mxy( ll, ss ), val.magNRMSEInPlaneCentralSlice_FA( ll, ss ),...
                val.meanFAInPlaneCentralSlice( ll, ss ), val.stdFAInPlaneCentralSlice( ll, ss ), val.nMaxMinInPlaneCentralSlice( ll, ss ), val.magCoeffVarInPlaneCentralSlice( ll, ss ) ] = ...
                calculateFAMtargMetrics( MvecInPlaneCentralSlice, FAvecInPlaneCentralSlice, MtargInPlaneCentralSlice, FAtargInPlaneCentralSlice );
        else
            val.calcInPlaneMetrics( ll, ss ) = false;
        end

    end

end

end