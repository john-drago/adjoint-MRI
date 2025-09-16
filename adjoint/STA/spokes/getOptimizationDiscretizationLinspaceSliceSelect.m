function [ ssdisc, sliceBounds ] = getOptimizationDiscretizationLinspaceSliceSelect(...
    withinSlice_di, outsideSlice_di, outsideSliceExt, sliceLocation, sliceThickness, sliceBounds, optbounds )

%% Initialize discretization for slice select direction
ssdisc = [];

sliceLocation = sort( sliceLocation );
numSlices = size( sliceLocation, 1 );

numZCoordInSlice = floor( sliceThickness / withinSlice_di );

for ss = 1:numSlices
     
     if ( ss == 1 ) && ( ss == numSlices )
        ss_less_st = optbounds( 1 );
        ss_less_end = sliceBounds( ss, 1 );
        ss_more_st = sliceBounds( ss, 2 );
        ss_more_end = optbounds( 2 );
     elseif ss == 1
        ss_less_st = optbounds( 1 );
        ss_less_end = sliceBounds( ss, 1 );
        ss_more_st = sliceBounds( ss, 2 );
        ss_more_end = mean( [ sliceBounds( ss, 2 ), sliceBounds( ss+1, 1 ) ] );
     elseif ss == numSlices
         ss_less_st = mean( [ sliceBounds( ss-1, 2 ), sliceBounds( ss, 1 ) ] );
         ss_less_end = sliceBounds( ss, 1 );
         ss_more_st = sliceBounds( ss, 2 );
         ss_more_end = optbounds( 2 );
     else
         ss_less_st = mean( [ sliceBounds( ss-1, 2 ), sliceBounds( ss, 1 ) ] );
         ss_less_end = sliceBounds( ss, 1 );
         ss_more_st = sliceBounds( ss, 2 );
         ss_more_end = mean( [ sliceBounds( ss, 2 ), sliceBounds( ss+1, 1 ) ] );
     end

     sliceSpacing = ( diff( sliceBounds( ss, : ) ) - ( numZCoordInSlice - 1 ) * withinSlice_di ) / 2 +...
         ( withinSlice_di * ( 0 : ( numZCoordInSlice- 1 ) ) );
     intraSliceDiscretization = sliceBounds( ss, 1 ) + sliceSpacing;

     outsideExtendLeft = fliplr( ( ss_less_end - outsideSliceExt * sliceThickness ) : -outsideSlice_di : ss_less_st );
     outsideSliceLeft = fliplr( ( ss_less_end ) : -withinSlice_di : ( ss_less_end - outsideSliceExt * sliceThickness ) );
     outsideSliceRight = ( ( ss_more_st ) : withinSlice_di : ( ss_more_st + outsideSliceExt * sliceThickness ) );
     outsideExtendRight = ( ( ss_more_st + outsideSliceExt * sliceThickness ) : outsideSlice_di : ss_more_end );

     ssdisc_ss = [...
         outsideExtendLeft,...
         outsideSliceLeft,...
         intraSliceDiscretization,...
         outsideSliceRight,...
         outsideExtendRight,...
         ];

     ssdisc = [ ssdisc, ssdisc_ss ]; %#ok
end

numDigAfterDecimal = 8;
ssdisc = round( ssdisc, numDigAfterDecimal );

ssdisc = unique( ssdisc );

end