function powerMetricFig = plotPowerMetrics( opt, binVis )

%% Initialize variables
fontSizeAx = 12;

%% Get Waveforms
wv = opt.generateWaveforms( opt.pOpt, opt );
wv.RF = complex( wv.breal, wv.bimag );
wv.RFphasor = complex( wv.brealphasor, wv.bimagphasor );

%% Get Power Metrices
% pulse power
pulseLength = wv.tvec( end ) + wv.dtvec( end ) / 2;
[ ~, ~, pulsePowerChannel ] = calcTotalMaxRFPower( ...
    wv.RFphasor, wv.dtvec, pulseLength, opt.Z0, opt.dutyCycle );
if opt.numXYCoils > 1
    plotPulsePower = true;
    maxPulsePowerChannel = max( pulsePowerChannel, [], 1 );
else
    plotPulsePower = false;
end

% local SAR
if isfield( opt, 'VOPs' )
    [ ~, ~, localSAR ] = calcLocalSAR( ...
        opt.VOPs, wv.RFphasor, wv.dtvec, pulseLength, opt.dutyCycle );
    plotLocalSAR = true;
    maxLocalSAR = max( localSAR, [], 1 );
else
    plotLocalSAR = false;
end

% global SAR
if isfield( opt, 'QGlobal' )
    [ ~, ~, globalSAR ] = calcGlobalSAR( ...
        opt.QGlobal, wv.RFphasor, wv.dtvec, pulseLength, opt.dutyCycle );
    plotGlobalSAR = true;
else
    plotGlobalSAR = false;
end

%% Determine how many plots
numPlots = 1;
if plotPulsePower
    numPlots = numPlots + 1;
end
if plotLocalSAR
    if opt.numXYCoils > 1
        numPlots = numPlots + 3;
    else
        numPlots = numPlots + 2;
    end
end
if plotGlobalSAR
    numPlots = numPlots + 1;
end


%% Define Params to Use During Plotting
lw = 2.5;
currPlot = 0;
ColorMax = 1;
ColorMin = 0.25;
OtherMax = 0.75;
rfColors = generateMultipleWaveformColormap( opt.numXYCoils, 'g', ColorMax, ColorMin, OtherMax );

%%  Plot the Param Output
figSize = [ 1 1 1100 1100 ];

if binVis
    powerMetricFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    powerMetricFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

PSDTL = tiledlayout( numPlots, 1,...
    'padding', 'compact');

%% RF power
[ tvec_box, pulsePowerChannel_box ] = makeBoxyWaveformAdjoint(...
    wv.tvec, wv.dtvec, pulsePowerChannel );

currPlot = currPlot + 1;
axRFpower = nexttile(PSDTL, currPlot, [1 1]);
pRFpower = plot( transpose( tvec_box )/1e-3, transpose( pulsePowerChannel_box ) );
set(pRFpower, 'linewidth', lw);
set(pRFpower, {'color'}, rfColors);
xlim( [0 tvec_box(end)/1e-3] );
grid(axRFpower, 'on');
axRFpower.TickLabelInterpreter='latex';
axRFpower.FontSize = fontSizeAx;
xticklabels( axRFpower, [] );
ylabel('Pulse Power (W)', 'interpreter', 'latex');
title(['\textbf{', 'Pulse Power per Channel', '}'], 'interpreter', 'latex');

if plotPulsePower
    [ tvec_box, maxPulsePowerChannel_box ] = makeBoxyWaveformAdjoint(...
        wv.tvec, wv.dtvec, maxPulsePowerChannel );

    currPlot = currPlot + 1;
    axRFmaxpower = nexttile(PSDTL, currPlot, [1 1]);
    pRFmaxpower = plot( transpose( tvec_box )/1e-3, transpose( maxPulsePowerChannel_box ) );
    set(pRFmaxpower, 'linewidth', lw);
    set(pRFmaxpower, 'color', "#0072BD");
    xlim( [0 tvec_box(end)/1e-3] );
    grid(axRFmaxpower, 'on');
    axRFmaxpower.TickLabelInterpreter='latex';
    axRFmaxpower.FontSize = fontSizeAx;
    xticklabels( axRFmaxpower, [] );
    ylabel('Pulse Power (W)', 'interpreter', 'latex');
    title(['\textbf{', 'Max Pulse Power', '}'], 'interpreter', 'latex');
end

%% Plot Local SAR
if plotLocalSAR
    currPlot = currPlot + 1;
    if opt.numXYCoils > 1
        axLocalSAR = nexttile(PSDTL, currPlot, [2 1]);
    else
        axLocalSAR = nexttile(PSDTL, currPlot, [1 1]);
    end
    imLSAR = imagesc( 'XData', wv.tvec/1e-3, 'YData', (1:size(opt.VOPs, 3)).',...
        'CData', localSAR);
    xlim( [0 tvec_box(end)/1e-3] );
    ylim( [0.5 size(opt.VOPs, 3) + 0.5] );
    axLocalSAR.TickLabelInterpreter='latex';
    axLocalSAR.FontSize = fontSizeAx;
    xticklabels( axLocalSAR, [] );
    yticks( axLocalSAR, axLocalSAR.YTick(2:(end-1)) )
    hotcmap = colorcet( 'L04' );
    colormap( axLocalSAR, hotcmap );
    ylabel('VOP Index', 'interpreter', 'latex');
    title(['\textbf{', 'Local SAR', '}'], 'interpreter', 'latex');

    cbimLSAR = colorbar;
    ylabel( cbimLSAR, 'Local SAR (W/kg)', 'interpreter', 'latex' )
    cbimLSAR.TickLabelInterpreter = 'latex';
    if opt.numXYCoils > 1
        currPlot = currPlot + 1;
    end

    [ tvec_box, maxLocalSAR_box ] = makeBoxyWaveformAdjoint(...
        wv.tvec, wv.dtvec, maxLocalSAR );

    currPlot = currPlot + 1;
    axMaxLocalSAR = nexttile(PSDTL, currPlot, [1 1]);
    pLocalSAR = plot( transpose( tvec_box )/1e-3, transpose( maxLocalSAR_box ) );
    set(pLocalSAR, 'linewidth', lw);
    set(pLocalSAR, 'color', 	"#A2142F");
    xlim( [0 tvec_box(end)/1e-3] );
    grid(axMaxLocalSAR, 'on');
    axMaxLocalSAR.TickLabelInterpreter='latex';
    axMaxLocalSAR.FontSize = fontSizeAx;
    xticklabels( axMaxLocalSAR, [] );
    ylabel('Max Local SAR (W/kg)', 'interpreter', 'latex');
    title(['\textbf{', 'Max Local SAR', '}'], 'interpreter', 'latex');

end

%% Plot Global SAR
if plotGlobalSAR
    [ tvec_box, globalSAR_box ] = makeBoxyWaveformAdjoint(...
        wv.tvec, wv.dtvec, globalSAR );

    currPlot = currPlot + 1;
    axGlobalSAR = nexttile(PSDTL, currPlot, [1 1]);
    pGSAR = plot( transpose( tvec_box )/1e-3, transpose( globalSAR_box ) );
    set(pGSAR, 'linewidth', lw);
    set(pGSAR, 'color', "#D95319" );
    xlim( [0 tvec_box(end)/1e-3] );
    grid(axGlobalSAR, 'on');
    axGlobalSAR.TickLabelInterpreter='latex';
    axGlobalSAR.FontSize = fontSizeAx;
    % xticklabels( axGlobalSAR, [] );
    ylabel('Global SAR (W/kg)', 'interpreter', 'latex');
    xlabel('Time (ms)', 'interpreter', 'latex')
    title(['\textbf{', 'Global SAR', '}'], 'interpreter', 'latex');

end

%% Set figure window style back to default
if ~binVis
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end