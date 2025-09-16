function [ opt, pulse ] = processAdjointTerminalCostFunction( opt, pulse )

if ~isempty( pulse.terminalCostFunction )
    if matches( pulse.terminalCostFunction, 'magnitude-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runMagnitudeLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'arcsin-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runArcsinLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'arccos-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runArccosLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'mxy-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runMxyLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'mxyz-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runMxyzLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'mz-least-squares', 'IgnoreCase', true  )
        opt.terminalCostFunction = @runMzLeastSquaresCost;
    elseif matches( pulse.terminalCostFunction, 'FA-with-mean-phase-entire-slice', 'IgnoreCase', true  )

        opt.terminalCostFunction = @runFAWithMeanPhaseEntireSliceCost;
        
        opt = getSliceIdxs( opt );

    elseif matches( pulse.terminalCostFunction, 'FA-with-mean-phase-across-slice', 'IgnoreCase', true  )

        opt.terminalCostFunction = @runFAWithMeanPhaseAcrossSliceCost;

        opt = getSliceIdxs( opt );

        [ opt.uniqueXYPos, opt.XYPosToUniqueXYPos, opt.uniqueXYPosToXYPos, opt.numUniqueXYPos  ] = ...
            getUniqueXYPos( opt.pos );
        
    else
        error( "Unknown 'pulse.terminalCostFunction'." );
    end
else
    opt.terminalCostFunction = [];
end

end