function [ oc, pulse, fields, optorval ] = processFieldsStruct(...
    oc, pulse, fields, si, optorval )
% This function will interpolate the fields to the discretization specified
% by di for the subject indices specified si and put them with the correct
% fieldname in the optorval struct.

%% Determine the correct grid size for optorval

if isstruct( optorval.di )

    if isfield( optorval.di, 'x' )
        optorval.x = optorval.di.x;
    end
    if isfield( optorval.di, 'y' )
        optorval.y = optorval.di.y;
    end
    if isfield( optorval.di, 'z' )
        optorval.z = optorval.di.z;
    end

elseif isnumeric( optorval.di )

    if isscalar( optorval.di )

        optorval.di = optorval.di * ones( 3, 1 );

    elseif numel( optorval.di ) == 3

        di = optorval.di( : );
        optorval.x = fields.x( 1 ) : di( 1 ) : fields.x( end );
        optorval.y = fields.y( 1 ) : di( 2 ) : fields.y( end );
        optorval.z = fields.z( 1 ) : di( 3 ) : fields.z( end );

    else
        error("Unsure how to process di for spacing.")
    end

else
    error("Unsure how to process di for spacing.")
end

if ( isfield( optorval, 'x' ) && isfield( optorval, 'y' ) ) && isfield( optorval, 'z' )
    [ optorval.X, optorval.Y, optorval.Z ] = ndgrid( optorval.x, optorval.y, optorval.z );
else
    error("Unsure how to process di for spacing.")
end

[ I, J, K ] = ndgrid( 1:length( optorval.x ), 1:length( optorval.y ), 1:length( optorval.z ) );

%% Determine number of coils
numZCoils = fields.numZCoils;
numXYCoils = fields.numXYCoils;
numSubj = length( si );

%% Initialize cell arrays to generate the needed variables
bzSens_vec_db = cell( numSubj, 1 );
b1p_vec_db = cell( numSubj, 1 );
db0_vec_db = cell( numSubj, 1 );
pos_vec_db = cell( numSubj, 1 );
idx_vec_db = cell( numSubj, 1 );
optroi_vec_db = cell( numSubj, 1 );
zi_vec_db = cell( numSubj, 1 );

M0_vec_db = cell( numSubj, 1 );
Mtarg_vec_db = cell( numSubj, 1 );

numPosSubj = zeros( numSubj, 1, 'uint32' );

subjIden = fields.subjIden( si );
subjPath = fields.subjPath( si );

%% Intepolate b1p
b1pInterpolate = zeros( [ numXYCoils, size( optorval.X ), numSubj ] );

if numSubj > 1
    for nn = 1:numXYCoils
        for ss = 1:numSubj
            b1pInterpolate( nn, :, :, :, ss ) = Interp3D(...
                fields.X, fields.Y, fields.Z, squeeze( fields.b1p( nn, :, :, :, si(ss) ) ),...
                optorval.X, optorval.Y, optorval.Z );
        end
    end
elseif numSubj == 1
    for nn = 1:numXYCoils
        b1pInterpolate( nn, :, :, : ) = Interp3D(...
            fields.X, fields.Y, fields.Z, squeeze( fields.b1p( nn, :, :, : ) ),...
            optorval.X, optorval.Y, optorval.Z );
    end
else
    error( "Unknown number of subjects." )
end

%% Interpolate db0
db0Interpolate = zeros( [ size( optorval.X ), numSubj ] );

for ss = 1:numSubj
    if numSubj > 1
        db0Interpolate( :, :, :, ss ) = Interp3D(...
            fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, :, si(ss) ) ),...
            optorval.X, optorval.Y, optorval.Z );
    elseif numSubj == 1
        db0Interpolate( :, :, : ) = Interp3D(...
            fields.X, fields.Y, fields.Z, squeeze( fields.db0( :, :, : ) ),...
            optorval.X, optorval.Y, optorval.Z );
    else
        error( "Unknown number of subjects." )
    end
end

%% Interpolate roi
optROIInterpolate = false( [ size( optorval.X ), numSubj ] );
valROIInterpolate = false( [ size( optorval.X ), numSubj ] );

maskInterpThresh = 0.50 - eps;

for ss = 1:numSubj

    if numSubj > 1
        optROIInterpolate( :, :, :, ss ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, :, si(ss) ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;

        valROIInterpolate( :, :, :, ss ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.val_roi( :, :, :, si(ss) ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;
    elseif numSubj == 1
        optROIInterpolate( :, :, : ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, : ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;

        valROIInterpolate( :, :, : ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.val_roi( :, :, : ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;
    else
        error( "Unknown number of subjects." )
    end
end

%% Process possible z types, either discrete coils or a surface
%% Process shim array
if pulse.shimarray
    % Interpolate bz
    if numZCoils > 0
        bzSensInterpolate = zeros( [ numZCoils, size( optorval.X ) ] );

        for nn = 1:numZCoils
            bzSensInterpolate( nn, :, :, : ) = Interp3D(...
                fields.X, fields.Y, fields.Z, squeeze( fields.bz( nn, :, :, : ) ),...
                optorval.X, optorval.Y, optorval.Z );
        end
    else
        bzSensInterpolate = zeros( [ 1, size( optorval.X ) ] );
    end
    
else
    error( "Unknown 'fields.bz' type." )
end

%% Interpolate M0
M0Interpolate = zeros( [ size( optorval.X ), 3 ] );

for dd = 1:3
    M0Interpolate( :, :, :, dd ) = Interp3D(...
        fields.X, fields.Y, fields.Z, squeeze( fields.M0( :, :, :, dd ) ),...
        optorval.X, optorval.Y, optorval.Z );
end

%% Interpolate Mtarg
MtargInterpolate = zeros( [ size( optorval.X ), 3 ] );

for dd = 1:3
    MtargInterpolate( :, :, :, dd ) = Interp3D(...
        fields.X, fields.Y, fields.Z, squeeze( fields.Mtarg( :, :, :, dd ) ),...
        optorval.X, optorval.Y, optorval.Z );
end

%% Fill in cell arrays
for ss = 1:numSubj

    % Define roi based on `opt` or `val`.
    switch optorval.structtype
        case {'opt'}
            roiSubj = optROIInterpolate( :, :, :, ss );
        case {'val'}
            roiSubj = valROIInterpolate( :, :, :, ss );
        otherwise
            error( "Unknown 'opt_or_val.structtype'." )
    end

    % Get index locations in 3D array for the roi.
    idx_vec_db{ ss } = uint32( find( roiSubj ) );

    % Find indices of opt roi (which is usually the brain roi)
    switch optorval.structtype
        case {'opt'}
            optroi_vec_db{ ss } = true( size( idx_vec_db{ ss } ) );
        case {'val'}
            roiOptSubj = optROIInterpolate( :, :, :, ss );
            optIdx = uint32( find( roiOptSubj ) );
            optroi_vec_db{ ss } = ismember( idx_vec_db{ ss }, optIdx );
    end

    % Get slices of array where there is opt roi
    zi_vec_db{ ss } = find( squeeze( sum( roiSubj, [ 1 2 ] ) ) );

    % Get bz in the roisubj
    if numZCoils > 0
        bzSensSubj = zeros( size( idx_vec_db{ ss }, 1 ), numZCoils );
        for cc = 1:numZCoils
            coilField = bzSensInterpolate( cc, :, :, : );
            bzSensSubj( :, cc ) = coilField( roiSubj );
        end
        bzSens_vec_db{ ss } = bzSensSubj;
    else
        bzSensSubj = zeros( size( idx_vec_db{ ss }, 1 ), 1 );
        bzSens_vec_db{ ss } = bzSensSubj;
    end

    % Get b1p in the roisubj
    b1pSubj = zeros( size( idx_vec_db{ ss }, 1 ), numXYCoils );
    for cc = 1:numXYCoils
        coilField = b1pInterpolate( cc, :, :, :, ss );
        b1pSubj( :, cc ) = coilField( roiSubj );
    end
    b1p_vec_db{ ss } = b1pSubj;

    % Get db0
    db0Subj = db0Interpolate( :, :, :, ss );
    db0_vec_db{ ss } = db0Subj( roiSubj );

    % Get locations
    pos_vec_db{ ss } = [...
        optorval.X( roiSubj ), optorval.Y( roiSubj ), optorval.Z( roiSubj ) ];

    % Get number of points
    numPosSubj( ss ) = size( pos_vec_db{ ss }, 1 );

    % Get initial magnetization
    M0_vec = zeros( numPosSubj( ss ), 3 );
    Mtarg_vec = zeros( numPosSubj( ss ), 3 );

    for dd = 1:3
        Mind = sub2ind( size( M0Interpolate ),...
            I( roiSubj ), J( roiSubj ), K( roiSubj ), dd * ones( numPosSubj( ss ), 1 ) );

        M0_vec( :, dd ) = M0Interpolate( Mind );
        Mtarg_vec( :, dd ) = MtargInterpolate( Mind );
    end
    M0_vec_db{ ss } = M0_vec;

    % Get Mtarg
    Mtarg_vec_db{ ss } = Mtarg_vec;

end

%% Determine if need to update db0 because of frequency offset

db0 = cell2mat( db0_vec_db( : ) );

if isfield( pulse, 'constantRotatingFrame' )

    optorval.constantRotatingFrame = pulse.constantRotatingFrame;

    if pulse.constantRotatingFrame

        % Define dwxy = w0 - wxy, where w0 is Larmor and wxy is applied pulse
        % frequency
        if isfield( pulse, 'mdfxy' )
            pulse.dfxy = - pulse.mdfxy;
            pulse.dwxy = pulse.dfxy * ( 2 * pi );
            pulse.mdwxy = pulse.mdfxy * ( 2 * pi );
        elseif isfield( pulse, 'mdwxy' )
            pulse.mdfxy = pulse.mdwxy / ( 2 * pi );
            pulse.dwxy = - pulse.mdwxy;
            pulse.dfxy = pulse.dwxy / ( 2 * pi );
        elseif isfield( pulse, 'dfxy' )
            pulse.dwxy = pulse.dfxy * ( 2 * pi );
            pulse.mdfxy = - pulse.dfxy;
            pulse.mdwxy = pulse.mdfxy * ( 2 * pi );
        elseif isfield( pulse, 'dwxy' )
            pulse.mdwxy = - pulse.dwxy;
            pulse.mdfxy = pulse.mdwxy / ( 2 * pi );
            pulse.dfxy = pulse.dwxy / ( 2 * pi );
        else
            pulse.mdwxy = 0;
            pulse.mdfxy = 0;
            pulse.dwxy = 0;
            pulse.dfxy = 0;
        end

        optorval.mdwxy = pulse.mdwxy;
        optorval.mdfxy = pulse.mdfxy;
        optorval.dwxy = pulse.dwxy;
        optorval.dfxy = pulse.dfxy;

    else % non-constant rotating frame

        if isfield( pulse, 'timing_dfxy' )
            optorval.timing_dfxy = pulse.timing_dfxy;
            pulse.timing_dwxy = pulse.timing_dfxy;
            optorval.timing_dwxy = pulse.timing_dwxy;

        elseif isfield( pulse, 'timing_dwxy' )
            optorval.timing_dwxy = pulse.timing_dwxy;
            % pulse.timing_dfxy = pulse.timing_dwxy;
            % optorval.timing_dfxy = pulse.timing_dfxy;
        else
            warning( "Did not specify pulse.timing_dwxy. Switching to constant frame at dwxy = 0." )
            pulse.mdwxy = 0;
            pulse.mdfxy = 0;
            pulse.dwxy = 0;
            pulse.dfxy = 0;

            optorval.mdwxy = pulse.mdwxy;
            optorval.mdfxy = pulse.mdfxy;
            optorval.dwxy = pulse.dwxy;
            optorval.dfxy = pulse.dfxy;

            pulse.constantRotatingFrame = true;
            optorval.constantRotatingFrame = true;
        end

        if isfield( pulse, 'vals_dfxy' )
            optorval.vals_dfxy = pulse.vals_dfxy( : );
            pulse.vals_dwxy = 2*pi * pulse.vals_dfxy( : );
            optorval.vals_dwxy = pulse.vals_dwxy;
        elseif isfield( pulse, 'vals_dwxy' )
            optorval.vals_dwxy = pulse.vals_dwxy( : );
        else
            warning( "Did not specify pulse.vals_dwxy. Switching to constant frame at dwxy = 0." )
            pulse.mdwxy = 0;
            pulse.mdfxy = 0;
            pulse.dwxy = 0;
            pulse.dfxy = 0;

            optorval.mdwxy = pulse.mdwxy;
            optorval.mdfxy = pulse.mdfxy;
            optorval.dwxy = pulse.dwxy;
            optorval.dfxy = pulse.dfxy;

            pulse.constantRotatingFrame = true;
            optorval.constantRotatingFrame = true;
        end

        if ~pulse.constantRotatingFrame
            if size( optorval.timing_dwxy, 1 ) == length( optorval.vals_dwxy )
            elseif size( optorval.timing_dwxy, 2 ) == length( optorval.vals_dwxy )
                optorval.timing_dwxy = transpose( optorval.timing_dwxy );
            else
                error( "Different number of dwxy periods and values." )
            end
        end

    end

else
    pulse.mdwxy = 0;
    pulse.mdfxy = 0;
    pulse.dwxy = 0;
    pulse.dfxy = 0;

    optorval.mdwxy = pulse.mdwxy;
    optorval.mdfxy = pulse.mdfxy;
    optorval.dwxy = pulse.dwxy;
    optorval.dfxy = pulse.dfxy;

    pulse.constantRotatingFrame = true;
    optorval.constantRotatingFrame = true;
end

%% Determine whether to change magnetization
if isfield( pulse, 'convertMBackToLarmor' )
    optorval.convertMBackToLarmor = pulse.convertMBackToLarmor;
else
    optorval.convertMBackToLarmor = false;
    pulse.convertMBackToLarmor = false;
end

if pulse.constantRotatingFrame && ( pulse.dwxy == 0 )
    optorval.convertMBackToLarmor = false;
    pulse.convertMBackToLarmor = false;
end

%% Determine whether to shift Mtarg
if isfield( pulse, 'convertMtargAwayLarmor' )
    optorval.convertMtargAwayLarmor = pulse.convertMtargAwayLarmor;
else
    optorval.convertMtargAwayLarmor = false;
    pulse.convertMtargAwayLarmor = false;
end

if pulse.constantRotatingFrame && ( pulse.dwxy == 0 )
    optorval.convertMtargAwayLarmor = false;
    pulse.convertMtargAwayLarmor = false;
end

%% Deal with VOPs and Global SAR if present
if isfield( fields, 'VOPs' )
    optorval.VOPs = fields.VOPs;
    optorval.numVOPs = size( optorval.VOPs, 3 );
end
if isfield( fields, 'QGlobal' )
    optorval.QGlobal = fields.QGlobal;
end
if isfield( fields, 'VOPs_E10m' )
    optorval.VOPs_E10m = fields.VOPs_E10m;
    optorval.numVOPs_E10m = size( optorval.VOPs_E10m, 3 );
end
if isfield( fields, 'VOPs_H10m' )
    optorval.VOPs_H10m = fields.VOPs_H10m;
    optorval.numVOPs_H10m = size( optorval.VOPs_H10m, 3 );
end

%% Add values to struct
optorval.bzsens = cell2mat( bzSens_vec_db( : ) );
optorval.b1p = cell2mat( b1p_vec_db( : ) );
optorval.b1preal = real( optorval.b1p );
optorval.b1pimag = imag( optorval.b1p );
optorval = rmfield( optorval, "b1p" );
optorval.db0 = db0;
optorval.pos = cell2mat( pos_vec_db( : ) );
optorval.M0 = cell2mat( M0_vec_db( : ) );
optorval.Mtarg = cell2mat( Mtarg_vec_db( : ) );
optorval.idx = cell2mat( idx_vec_db( : ) );
optorval.optroi = cell2mat( optroi_vec_db( : ) );

optorval.zi = zi_vec_db; % cell array of the indices where roi is

gyro = 267.5e6;
optorval.gyro = gyro;

optorval.numPos = size( optorval.M0, 1 );
optorval.numPosSubj = numPosSubj;
optorval.cum_numPosSubj = cumsum( optorval.numPosSubj );

optorval.subjIden = subjIden;
optorval.subjPath = subjPath;

optorval.numZCoils = numZCoils;
optorval.numXYCoils = numXYCoils;
optorval.numSubj = numSubj;

end