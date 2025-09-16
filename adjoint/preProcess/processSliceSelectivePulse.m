function [ opt, pulse ] = processSliceSelectivePulse( opt, pulse )

%% Sort the slice Bounds
sliceBounds = zeros( size( pulse.sliceBounds ) );
[ sliceBounds( :, 1 ), leftBoundsSortIdx ] = sort( pulse.sliceBounds( :, 1 ), 'ascend' );
[ sliceBounds( :, 2 ), rightBoundsSortIdx ] = sort( pulse.sliceBounds( :, 2 ), 'ascend' );

[ ~, sliceBoundsMeanIdx ] = sort( mean( sliceBounds, 2 ), 'ascend' );

[ sliceLocation, sliceLocationIdx ] = sort( pulse.sliceLocation, 'ascend' );

if ~all( sliceLocationIdx == sliceBoundsMeanIdx )
    error( "Error with the ordering of locations and bounds" )
end

if ~all( leftBoundsSortIdx == rightBoundsSortIdx )
    error( "Slices overlap with each other" )
end

%% Determine slice indices from
% opt.sliceIdxs = zeros( pulse.numSlices, 2  );
opt.slicesMIdx = false( opt.numPos, pulse.numSlices );

if isfield( opt, 'z' )
    ditol = min( abs( diff( opt.z ) ) ) / 10;
else
    ditol = 5e-5;
end

for ss = 1:pulse.numSlices
    % opt.sliceIdxs( ss, 1 ) = find( ( opt.z >= ( sliceBounds( ss, 1 ) - ditol ) ), 1, 'first' );
    % opt.sliceIdxs( ss, 2 ) = find( ( opt.z <= ( sliceBounds( ss, 2 ) + ditol ) ), 1, 'last' );

    opt.slicesMIdx( :, ss ) = ...
        ( opt.pos( :, 3 ) >= ( sliceBounds( ss, 1 ) - ditol ) ) &...
        ( opt.pos( :, 3 ) <= ( sliceBounds( ss, 2 ) + ditol ) );
end

opt.allSlicesMIdx = any( opt.slicesMIdx, 2 );

%% Assign to pulse struct
pulse.sliceBounds = sliceBounds;
pulse.sliceLocation = sliceLocation;

%% Assign to opt struct
opt.numSlices = pulse.numSlices;
opt.sliceLocation = sliceLocation;
opt.sliceBounds = sliceBounds;
opt.sliceThickness = pulse.sliceThickness;

end