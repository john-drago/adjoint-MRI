function [ cellFigs, strLabels ] = processAdjointImages( opt, pulse, valtrain, valtest, oc )

%% Initialize cell array
cellFigs = cell( 1, 1 );
strLabels = strings( 1, 1 );
imageidx = 1;

%% Get opt train figures
[ optTrainFAFig, optTrainPhFig ] = plotResultAdjoint(...
        valtrain, pulse, "opt", oc.binVis );

cellFigs{ imageidx, 1 } = optTrainFAFig;
strLabels( imageidx, 1 ) = "optTrainFAFig";

imageidx = imageidx + 1;

cellFigs{ imageidx, 1 } = optTrainPhFig;
strLabels( imageidx, 1 ) = "optTrainPhFig";

imageidx = imageidx + 1;

%% Get opt test figures
if ~isempty( valtest )
    [ optTestFAFig, optTestPhFig ] = plotResultAdjoint(...
        valtest, pulse, "opt", oc.binVis );

    cellFigs{ imageidx, 1 } = optTestFAFig;
    strLabels( imageidx, 1 ) = "optTestFAFig";
    imageidx = imageidx + 1;

    cellFigs{ imageidx, 1 } = optTestPhFig;
    strLabels( imageidx, 1 ) = "optTestPhFig";
    imageidx = imageidx + 1;
end

%% Generate Slice Images
if isfield( pulse, 'sliceSelective' ) && pulse.sliceSelective

    [ thruSliceTrainFAFig, thruSliceTrainPhFig ] = plotThruSlice( valtrain, oc.binVis );

    cellFigs{ imageidx, 1 } = thruSliceTrainFAFig;
    strLabels( imageidx, 1 ) = "thruSliceTrainFA";
    imageidx = imageidx + 1;

    cellFigs{ imageidx, 1 } = thruSliceTrainPhFig;
    strLabels( imageidx, 1 ) = "thruSliceTrainPh";
    imageidx = imageidx + 1;

    % [ inSliceFATrain, inSlicePhTrain ] = plotInSlice( valtrain, oc.binVis );
    % 
    % cellFigs{ imageidx, 1 } = inSliceFATrain;
    % strLabels( imageidx, 1 ) = "inSliceFATrain";
    % imageidx = imageidx + 1;
    % 
    % cellFigs{ imageidx, 1 } = inSlicePhTrain;
    % strLabels( imageidx, 1 ) = "inSlicePhTrain";
    % imageidx = imageidx + 1;

    if ~isempty( valtest )
        
        [ thruSliceTestFAFig, thruSliceTestPhFig ] = plotThruSlice( valtest, oc.binVis );

        cellFigs{ imageidx, 1 } = thruSliceTestFAFig;
        strLabels( imageidx, 1 ) = "thruSliceTestFA";
        imageidx = imageidx + 1;

        cellFigs{ imageidx, 1 } = thruSliceTestPhFig;
        strLabels( imageidx, 1 ) = "thruSliceTestPhFig";
        imageidx = imageidx + 1;

        % [ inSliceFATest, inSlicePhTest ] = plotInSlice( valtest, oc.binVis );
        % 
        % cellFigs{ imageidx, 1 } = inSliceFATest;
        % strLabels( imageidx, 1 ) = "inSliceFATest";
        % imageidx = imageidx + 1;
        % 
        % cellFigs{ imageidx, 1 } = inSlicePhTest;
        % strLabels( imageidx, 1 ) = "inSlicePhTest";
        % imageidx = imageidx + 1;
    end
end

%% Get Birdcage figures
if ( valtrain.numXYCoils == 1 ) && ( valtrain.performBCHPcomp )

    [ birdcageTrainFAFig, birdcageTrainPhFig ] = plotResultAdjoint(...
        valtrain, pulse, "BCHP", oc.binVis );

    cellFigs{ imageidx, 1 } = birdcageTrainFAFig;
    strLabels( imageidx, 1 ) = "birdcageTrainFAFig";
    imageidx = imageidx + 1;

    cellFigs{ imageidx, 1 } = birdcageTrainPhFig;
    strLabels( imageidx, 1 ) = "birdcageTrainPhFig";
    imageidx = imageidx + 1;

    if ~isempty( valtest )

        [ birdcageTestFAFig, birdcageTestPhFig ] = plotResultAdjoint(...
            valtest, pulse, "BCHP", oc.binVis );

        cellFigs{ imageidx, 1 } = birdcageTestFAFig;
        strLabels( imageidx, 1 ) = "birdcageTestFAFig";
        imageidx = imageidx + 1;

        cellFigs{ imageidx, 1 } = birdcageTestPhFig;
        strLabels( imageidx, 1 ) = "birdcageTestPhFig";
        imageidx = imageidx + 1;
    end
end

%% Generate Pulse Sequence Diagram
cellFigs{ imageidx, 1 } = plotAdjointPSD( valtrain, oc.binVis );
strLabels( imageidx, 1 ) = "psdFig";

imageidx = imageidx + 1;

%% Generate Pulse Sequence Diagram of the optimization (boxy waveform)
cellFigs{ imageidx, 1 } = plotOptWaveforms( opt, oc.binVis );
strLabels( imageidx, 1 ) = "optWaveformFig";

imageidx = imageidx + 1;

%% Generate Diagram of Slew Rate and Accel for RF, Grad, Shim
cellFigs{ imageidx, 1 } = plotSlewAccelWaveforms( opt, oc.binVis );
strLabels( imageidx, 1 ) = "slewAccelWaveformFig";

imageidx = imageidx + 1;

%% Generate Diagram of Power Metrics
cellFigs{ imageidx, 1 } = plotPowerMetrics( opt, oc.binVis );
strLabels( imageidx, 1 ) = "powerMetricFig";

imageidx = imageidx + 1;

%% Generate Freq Figures
if ( isfield( pulse, 'sliceSelective' ) && pulse.sliceSelective ) ||...
        any( strcmpi( pulse.name, [ "pp"; "cheb"; "pwc"; "fourier" ] ) )
    cellFigs{ imageidx, 1 } = plotPulseFreq( valtrain, oc.binVis );
    strLabels( imageidx, 1 ) = "RFFreq";

    imageidx = imageidx + 1;
end

%% Plot Convergence Plots
if opt.trackConvergence
    
    convFig = plotOptConvergence( opt, oc.binVis );

    cellFigs{ imageidx, 1 } = convFig;
    strLabels( imageidx, 1 ) = "convergenceFig";
    imageidx = imageidx + 1; %#ok

end

%% Save Images
if oc.saveResult
    imgSavePath = fullfile( oc.saveDir, "figs" );
    saveTypeFlag = 2; % 1 is .png only, 2 is with .fig, 3 is with .eps
    saveMultipleImagesAdjoint( imgSavePath, strLabels, cellFigs, oc.binVis, saveTypeFlag );
end

end