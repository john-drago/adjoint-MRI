function [ oc, pulse, fields ] = processFieldsData( oc, pulse, fields  )
% This function will process the fields struct to understand the geometry
% of the problem

%% Determine pulse type
switch pulse.optPopulation
    case {'tailored'}
        numSubj = 1;
        pulse.universal = false;
        pulse.tailored = true;
        if length( fields.subjIden ) > 1
            error( "Too many subjects for tailored pulse design." )
        end
    case {'universal'}
        numSubj = length( fields.subjIden );
        pulse.universal = true;
        pulse.tailored = false;
    otherwise
        error( "Unknown 'pulse.optPopulation'." )
end

%% Process indices of the subject population if tailored or universal
switch pulse.optPopulation
    case {'tailored'}
        si_train = 1; % subject index train
        fields.si_train = si_train;
        fields.si_test = [];
    case {'universal'}
        if isa( pulse.universalTrain, 'numeric' ) &&...
                isa( pulse.universalTest, 'numeric' )
            fields.si_train = uint32( pulse.universalTrain );
            fields.si_test = uint32( pulse.universalTest );
        elseif isa( pulse.universalTrain, 'logical' ) &&...
                isa( pulse.universalTest, 'logical' )
            fields.si_train = uint32( find( pulse.universalTrain ) );
            fields.si_test = uint32( find( pulse.universalTest ) );
        elseif isa( pulse.universalTrain, 'string' ) &&...
                isa( pulse.universalTest, 'string' )
            fields.si_train = uint32( find( contains( fields.subjIden, pulse.universalTrain, 'ignorecase', true ) ) );
            fields.si_test = uint32( find( contains( fields.subjIden, pulse.universalTest, 'ignorecase', true ) ) );
        elseif isa( pulse.universalTrain, 'cell' ) &&...
                isa( pulse.universalTest, 'cell' )
            pulse.universalTrain = string( pulse.universalTrain );
            pulse.universalTest = string( pulse.universalTest );
            fields.si_train = uint32( find( contains( fields.subjIden, pulse.universalTrain, 'ignorecase', true ) ) );
            fields.si_test = uint32( find( contains( fields.subjIden, pulse.universalTest, 'ignorecase', true ) ) );
        else
            error( "Unknown" )
        end
    otherwise
        error( "Unknown 'pulse.optPopulation'." )
end

%% Determine type of z-directed coil (if discrete windings) or if surface
if isnumeric( fields.bz )

    pulse.shimarray = true;

    oc.shimarray = true;

    bzSize = size( fields.bz );

    if ( isscalar( fields.bz ) ) && ( fields.bz == 0 )
        numZCoils = 0;
    elseif length( bzSize ) == 3
        numZCoils = 1;
        fields.bz = reshape( fields.bz, [ 1, bzSize ] );
    elseif length( bzSize ) == 4
        numZCoils = bzSize( 1 );
    else
        error( "Can't process size of 'fields.bz'." )
    end

else
    error( "Unknown 'fields.bz' type." )
end

%% Determine number of xy-coils
bxySize = size( fields.b1p );

if length( bxySize ) == 3
    numXYCoils = 1;
    fields.b1p = reshape( fields.b1p, [ 1, bxySize, 1 ] );
elseif length( bxySize ) == 4
    if bxySize(end) == numSubj % universal single channel
        numXYCoils = 1;
        fields.b1p = reshape( fields.b1p, [ 1, bxySize(1:3), numSubj ] );
    elseif all( bxySize( 2:end ) == size( fields.X ) ) % tailored pTx
        numXYCoils = bxySize( 1 );
    else
        error( "Can't process size of 'fields.b1p'." )
    end
elseif ( length( bxySize ) == 5 ) && ( numSubj == bxySize( end ) )
    numXYCoils = bxySize( 1 );
else
    error( "Can't process size of 'fields.b1p'." )
end

%% Refine size of db0
db0Size = size( fields.db0 );

if length( db0Size ) == 3
    % fields.db0 = reshape( fields.db0, [ 1, db0Size ] ); % db0 is correct size
elseif ( length( db0Size ) == 4 ) && ( db0Size(end) == numSubj )
else
    error( "Can't process size of 'fields.db0'." )
end

%% Refine size of opt_roi
opt_roi_Size = size( fields.opt_roi );

if length( opt_roi_Size ) == 3
    % fields.opt_roi = reshape( fields.opt_roi, [ 1, opt_roi_Size ] );
elseif ( length( opt_roi_Size ) == 4 ) && ( opt_roi_Size(end) == numSubj )
else
    error( "Can't process size of 'fields.opt_roi'." )
end

%% Refine size of val_roi
val_roi_Size = size( fields.val_roi );

if length( val_roi_Size ) == 3
    % fields.val_roi = reshape( fields.val_roi, [ 1, val_roi_Size ] );
elseif ( length( val_roi_Size ) == 4 ) && ( val_roi_Size(end) == numSubj )
else
    error( "Can't process size of 'fields.val_roi'." )
end

%% Determine MTarg Array
[ fields.Mtarg, fields.targSt ] = generateMtargArray( pulse.targSt, pulse.targSt.MtargVals );

%% Determine M0 Array
[ fields.M0, fields.targSt ] = generateM0Array( pulse.targSt, pulse.targSt.M0Vals );

%% Determine spacing for Mtarg and M0 to see if need to interpolate the fields 
if ( numel( pulse.targSt.x ) == numel( fields.x ) ) && ...
        ( numel( pulse.targSt.y ) == numel( fields.y ) ) &&...
        ( numel( pulse.targSt.z ) == numel( fields.z ) )
    dEps = eps( 1e3 );
    allx = all( abs( pulse.targSt.x - fields.x ) < dEps );
    ally = all( abs( pulse.targSt.y - fields.y ) < dEps );
    allz = all( abs( pulse.targSt.z - fields.z ) < dEps );

    alld = ( allx & ally ) & allz;
else
    alld = false;
end

%% Interpolate the fields at the resolution specified by the MTarg and M0 array
if ~alld

    domSize = size( fields.Mtarg, [ 1, 2, 3 ] );

    % First b1p
    b1pInterpolate = zeros( [ numXYCoils, domSize, numSubj ] );
    if numSubj > 1
        for nn = 1:numXYCoils
            for ss = 1:numSubj
                b1pInterpolate( nn, :, :, :, ss ) = Interp3D(...
                    fields.X, fields.Y, fields.Z, squeeze( fields.b1p( nn, :, :, :, ss ) ),...
                    fields.targSt.X, fields.targSt.Y, fields.targSt.Z );
            end
        end
    elseif numSubj == 1
        for nn = 1:numXYCoils
            b1pInterpolate( nn, :, :, : ) = Interp3D(...
                fields.X, fields.Y, fields.Z, squeeze( fields.b1p( nn, :, :, : ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z );
        end
    else
        error( "Unknown number of subjects." )
    end

    % db0
    db0Interpolate = zeros( [ domSize, numSubj ] );

    for ss = 1:numSubj
        if numSubj > 1
            db0Interpolate( :, :, :, ss ) = Interp3D(...
                fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, :, ss ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z );
        elseif numSubj == 1
            db0Interpolate( :, :, : ) = Interp3D(...
                fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, : ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z );
        else
            error( "Unknown number of subjects." )
        end
    end

    % roi
    optROIInterpolate = false( [ domSize, numSubj ] );
    valROIInterpolate = false( [ domSize, numSubj ] );

    maskInterpThresh = 0.50 - eps;

    for ss = 1:numSubj

        if numSubj > 1
            optROIInterpolate( :, :, :, ss ) = Interp3D(...
                fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, :, ss ) ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z ) > maskInterpThresh;

            valROIInterpolate( :, :, :, ss ) = Interp3D(...
                fields.X, fields.Y, fields.Z, double( squeeze( fields.val_roi( :, :, :, ss ) ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z ) > maskInterpThresh;
        elseif numSubj == 1
            optROIInterpolate( :, :, : ) = Interp3D(...
                fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, : ) ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z ) > maskInterpThresh;

            valROIInterpolate( :, :, : ) = Interp3D(...
                fields.X, fields.Y, fields.Z, double( squeeze( fields.val_roi( :, :, : ) ) ),...
                fields.targSt.X, fields.targSt.Y, fields.targSt.Z ) > maskInterpThresh;
        else
            error( "Unknown number of subjects." )
        end
    end

    % process possible z types
    % shim array
    if pulse.shimarray
        % Interpolate bz
        if numZCoils > 0
            bzSensInterpolate = zeros( [ numZCoils, domSize ] );

            for nn = 1:numZCoils
                bzSensInterpolate( nn, :, :, : ) = Interp3D(...
                    fields.X, fields.Y, fields.Z, squeeze( fields.bz( nn, :, :, : ) ),...
                    fields.targSt.X, fields.targSt.Y, fields.targSt.Z );
            end
        else
            bzSensInterpolate = zeros( [ 1, domSize ] );
        end

    else
        error( "Unknown 'fields.bz' type." )
    end
end

%% Put new fields into the fields struct
if ~alld
    fields.X = fields.targSt.X;
    fields.Y = fields.targSt.Y;
    fields.Z = fields.targSt.Z;
    fields.b1p = b1pInterpolate;
    fields.db0 = db0Interpolate;
    fields.opt_roi = optROIInterpolate;
    fields.val_roi = valROIInterpolate;
    if pulse.shimarray
        if numZCoils > 0
            fields.bz = bzSensInterpolate;
        else
            fields.bz = 0;
        end
    end
end

%% Add parameters to structs
fields.numXYCoils = numXYCoils;
fields.numZCoils = numZCoils;

%% Get info about fields to pulse struct
pulse.subjectMaps = fields.subjectMaps;
if isfield( fields, 'zCoilPath' )
    pulse.zCoilPath = fields.zCoilPath;
else
    pulse.zCoilPath = [];
end
if isfield( fields, 'zCoilShiftPath' )
    pulse.zCoilShiftPath = fields.zCoilShiftPath;
else
    pulse.zCoilShiftPath = [];
end
if isfield( fields, 'VOPpath' )
    pulse.VOPpath = fields.VOPpath;
else
    pulse.VOPpath = [];
end

pulse.subjIden = fields.subjIden;
pulse.subjPath = fields.subjPath;

end