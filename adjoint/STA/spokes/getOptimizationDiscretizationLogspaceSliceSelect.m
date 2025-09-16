function [ ssdisc, sliceBounds ] = getOptimizationDiscretizationLogspaceSliceSelect(...
    withinSlice_di, outsideSlice_di_start, outsideSlice_di_end,...
    outsideSliceExt, sliceLocation, sliceThickness, sliceBounds, optbounds )

%% Initialize discretization for slice select direction
ssdisc = [];

sliceLocation = sort( sliceLocation );
numSlices = size( sliceLocation, 1 );

numZCoordInSlice = floor( sliceThickness / withinSlice_di );

baseMultFactor = 1.05;
minPoints = 20;
outsideSlice_di_steps = [...
    10.^( log10( outsideSlice_di_start ) : log10( baseMultFactor ) : log10( outsideSlice_di_end ) ),...
    outsideSlice_di_end * ones( 1, ceil( abs( diff( optbounds ) )/outsideSlice_di_end ) ) ];
outsideSlice_di_cumsum = cumsum( outsideSlice_di_steps );

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

     outsideExtendLeftAmt = abs( ( ss_less_end - outsideSliceExt * sliceThickness ) - ss_less_st );
     outsideExtendLeft = fliplr( ( ss_less_end - outsideSliceExt * sliceThickness ) ...
         - outsideSlice_di_cumsum( outsideSlice_di_cumsum < outsideExtendLeftAmt ) );
     if ( length( outsideExtendLeft ) < minPoints )
         outsideExtendLeft = fliplr( ( ss_less_end - outsideSliceExt * sliceThickness ) ...
         - logspace( log10(outsideSlice_di_start), log10( outsideExtendLeftAmt ), minPoints ) );
     end

     outsideSliceLeft = fliplr( ( ss_less_end ) : -withinSlice_di : ( ss_less_end - outsideSliceExt * sliceThickness ) );
     outsideSliceRight = ( ( ss_more_st ) : withinSlice_di : ( ss_more_st + outsideSliceExt * sliceThickness ) );

     outsideExtendRightAmt = abs( ss_more_end - ( ss_more_st + outsideSliceExt * sliceThickness ) );
     outsideExtendRight = ( ss_more_st + outsideSliceExt * sliceThickness ) ...
         + outsideSlice_di_cumsum( outsideSlice_di_cumsum < outsideExtendRightAmt );

     if ( length( outsideExtendLeft ) < minPoints )
         outsideExtendRight = ( ss_more_st + outsideSliceExt * sliceThickness ) ...
             + logspace( log10(outsideSlice_di_start), log10( outsideExtendRightAmt ), minPoints );
     end

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