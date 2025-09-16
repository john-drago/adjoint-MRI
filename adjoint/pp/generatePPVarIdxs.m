function [ PPVarIdxs, PPSFIdxs, PPVarSFIdxs, PPVarIdxsCell, PPVarIdxsCoilCell, numVarsPerChannel ] =...
    generatePPVarIdxs( numPPPeriods, numPPShape, numCoils, PPStartStop, idxtol )

if nargin < 5
    idxtol = 5e-7;
end

PPVarIdxs = zeros( numPPPeriods, 2, "uint32" );
PPSFIdxs = zeros( numPPPeriods, 2, "uint32" );
PPVarSFIdxs = zeros( numPPPeriods, 2, "uint32" );
PPVarIdxsCell = cell( numPPPeriods, 1 );
PPVarIdxsCoilCell = cell( numPPPeriods, 1 );

numVarsPerChannel = 0;

for pp = 1:numPPPeriods

    if ( pp == 1 ) && ( pp == numPPPeriods )
        PPVarIdxs( pp, : ) = [ 1, (numPPShape-2) ];
        PPVarSFIdxs( pp, : ) = [ 2, numPPShape-1 ];
        numVarsPerChannel = numVarsPerChannel + (numPPShape-2);
    elseif pp == 1
        zeroSFRightBound = determineZeroShapeFunctionRightBoundary( pp, PPStartStop, idxtol );
        if zeroSFRightBound
            numVarsPerChannel = numVarsPerChannel + (numPPShape-1);
            PPVarIdxs( pp, : ) = [ 1, (numPPShape-1) ];
            PPVarSFIdxs( pp, : ) = [ 2, numPPShape ];
        else
            numVarsPerChannel = numVarsPerChannel + (numPPShape-2);
            PPVarIdxs( pp, : ) = [ 1, (numPPShape-2) ];
            PPVarSFIdxs( pp, : ) = [ 2, numPPShape-1 ];
        end
    elseif pp < numPPPeriods
        zeroSFRightBound = determineZeroShapeFunctionRightBoundary( pp, PPStartStop, idxtol );
        if prevZeroRightBound && zeroSFRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 0, numPPShape-1 ];
            PPVarSFIdxs( pp, : ) = [ 1, numPPShape ];
        elseif ~prevZeroRightBound && ~zeroSFRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 1, numPPShape-2 ];
            PPVarSFIdxs( pp, : ) = [ 2, numPPShape-1 ];
        elseif ~prevZeroRightBound && zeroSFRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 1, numPPShape-1 ];
            PPVarSFIdxs( pp, : ) = [ 2, numPPShape ];
        elseif prevZeroRightBound && ~zeroSFRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 0, numPPShape-2 ];
            PPVarSFIdxs( pp, : ) = [ 1, numPPShape-1 ];
        end
        if zeroSFRightBound
            numVarsPerChannel = numVarsPerChannel + (numPPShape-1);
        else
            numVarsPerChannel = numVarsPerChannel + (numPPShape-2);
        end

    elseif pp == numPPPeriods

        if prevZeroRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 0, numPPShape-2 ];
            PPVarSFIdxs( pp, : ) = [ 1, numPPShape-1 ];
        elseif ~prevZeroRightBound
            PPVarIdxs( pp, : ) = PPVarIdxs( pp-1, 2 ) +...
                [ 1, numPPShape-2 ];
            PPVarSFIdxs( pp, : ) = [ 2, numPPShape-1 ];
        end
        numVarsPerChannel = numVarsPerChannel + (numPPShape-2);
    end

    PPSFIdxs( pp, : ) = ( pp-1 ) * numPPShape + [ 1, numPPShape ];
    
    PPVarIdxsCell{ pp } = transpose( PPVarIdxs( pp, 1 ):PPVarIdxs( pp, 2 ) );

    prevZeroRightBound = zeroSFRightBound;
    
end

if numVarsPerChannel ~= PPVarIdxs( end, 2 )
    error( "Something went wrong with the variable assignment process." )
end

for pp = 1:numPPPeriods

    PPVarIdxsCoilCell{ pp } = repmat( PPVarIdxsCell{ pp },...
        [ 1, numCoils] ) +...
            ((0:uint32(numCoils-1))*(numVarsPerChannel));
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function zeroSFRightBound = determineZeroShapeFunctionRightBoundary( pp, PPStartStop, idxtol )
zeroSFRightBound = abs( PPStartStop( pp, 2 ) - PPStartStop( pp+1, 1 ) ) < idxtol;
end
% ----------------------------------------------------------------------- %