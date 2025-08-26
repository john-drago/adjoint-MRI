function [ Mtargvec, MtargInterpolate ] = generateMtargOpt( MtargArr, si, fields, optorval )

MtargSize = size( MtargArr );

if isscalar( MtargArr )
    Mtarg = zeros( [ size( fields.X ), 3 ] );
    Mtarg( :, :, :, 1 ) = MtargArr;
elseif numel( MtargArr ) == 3
    Mtarg = zeros( [ size( fields.X ), 3 ] );
    Mtarg( :, :, :, 1 ) = MtargArr( 1 );
    Mtarg( :, :, :, 2 ) = MtargArr( 2 );
    Mtarg( :, :, :, 3 ) = MtargArr( 3 );
elseif ( length( MtargSize ) == 3 ) && all( MtargSize(1:3) == size( fields.X ) )
    Mtarg = zeros( [ size( fields.X ), 3 ] );
    Mtarg( :, :, :, 1 ) = MtargArr;
elseif length( MtargSize ) == 4 && all( MtargSize(1:3) == size( fields.X ) )
    Mtarg = MtargArr;
else
    error( "Can't process size of 'pulse.Mtarg'." )
end

numPosSubj = optorval.numPosSubj;
numSubj = optorval.numSubj;

Mtarg_vec_db = cell( numSubj, 1 );

[ I, J, K ] = ndgrid( 1:length( optorval.x ), 1:length( optorval.y ), 1:length( optorval.z ) );

optROIInterpolate = false( [ size( optorval.X ), numSubj ] );

maskInterpThresh = 0.50 - eps;

for ss = 1:numSubj

    if numSubj > 1
        optROIInterpolate( :, :, :, ss ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, :, si(ss) ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;
    elseif numSubj == 1
        optROIInterpolate( :, :, : ) = Interp3D(...
            fields.X, fields.Y, fields.Z, double( squeeze( fields.opt_roi( :, :, : ) ) ),...
            optorval.X, optorval.Y, optorval.Z ) > maskInterpThresh;
    else
        error( "Unknown number of subjects." )
    end

end

MtargInterpolate = zeros( [ size( optorval.X ), 3 ] );

for dd = 1:3
    MtargInterpolate( :, :, :, dd ) = Interp3D(...
        fields.X, fields.Y, fields.Z, squeeze( Mtarg( :, :, :, dd ) ),...
        optorval.X, optorval.Y, optorval.Z );
end

for ss = 1:numSubj

    roiSubj = optROIInterpolate( :, :, :, ss );
    
    Mtargvec_iter = zeros( numPosSubj( ss ), 3 );

    for dd = 1:3
        Mind = sub2ind( size( MtargInterpolate ),...
            I( roiSubj ), J( roiSubj ), K( roiSubj ), dd * ones( numPosSubj( ss ), 1 ) );
        
        Mtargvec_iter( :, dd ) = MtargInterpolate( Mind );
    end

    % Get Mtarg
    Mtarg_vec_db{ ss } = Mtargvec_iter;
end



Mtargvec = cell2mat( Mtarg_vec_db( : ) );

end