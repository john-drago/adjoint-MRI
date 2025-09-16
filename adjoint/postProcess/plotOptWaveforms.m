function psdFig = plotOptWaveforms( opt, binVis )

%% Initialize variables
fontSizeAx = 12;

%% Get Waveforms
wv = opt.generateWaveforms( opt.pOpt, opt );

%% Define Params to Use During Plotting

numZCoils = opt.numZCoils;
numXYCoils = opt.numXYCoils;

lw = 2.5;

ColorMax = 1;
ColorMin = 0.25;
OtherMax = 0.75;

shimColors = generateMultipleWaveformColormap( numZCoils, 'b', ColorMax, ColorMin, OtherMax );
rfColorsReal = generateMultipleWaveformColormap( numXYCoils, 'g', ColorMax, ColorMin, OtherMax );
rfColorsImag = generateMultipleWaveformColormap( numXYCoils, 'r', ColorMax, ColorMin, OtherMax );

gradColors = num2cell(lines(4), 2);
gradColors = gradColors(2:end);

%%  Plot the Param Output
figSize = [ 1 1 1100 500 ];

if binVis
    psdFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    psdFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

PSDTL = tiledlayout( 4, 1,...
    'padding', 'compact');

%% brealphasor plot
[ tvec, brealphasor ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.breal );
% [ tvec, brealphasor ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.brealphasor );

axRFbreal = nexttile(PSDTL, 1);
pRF = plot( transpose( tvec )/1e-3, transpose( brealphasor ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsReal);
xlim( [0 tvec(end)/1e-3] );
ylim( [min(brealphasor, [], 'all'), max(brealphasor, [], 'all')] );
grid(axRFbreal, 'on');
axRFbreal.TickLabelInterpreter='latex';
axRFbreal.FontSize = fontSizeAx;
xticklabels( axRFbreal, [] );
ylabel('$b_{\rm real}$ (V)', 'interpreter', 'latex');
title(['\textbf{', 'RF real', '}'], 'interpreter', 'latex');

%% bimagphasor plot
[ tvec, bimagphasor ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.bimag );
% [ tvec, bimagphasor ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.bimagphasor );

axRFimag = nexttile(PSDTL, 2);
pRF = plot( transpose( tvec )/1e-3, transpose( bimagphasor ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsImag);
xlim( [0 tvec(end)/1e-3] );
ylim( [min(bimagphasor, [], 'all'), max(bimagphasor, [], 'all')] );
grid(axRFimag, 'on');
axRFimag.TickLabelInterpreter='latex';
axRFimag.FontSize = fontSizeAx;
xticklabels( axRFimag, [] );
ylabel('$b_{\rm imag}$ (V)', 'interpreter', 'latex');
title(['\textbf{', 'RF imag', '}'], 'interpreter', 'latex');

%% Shim Pulse Diagram
[ tvec, Shim ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.Shim );

axShimwaveform = nexttile(PSDTL, 3);
pShim = plot( transpose( tvec )/1e-3,  transpose( Shim ) );
set(pShim, 'linewidth', lw);
set(pShim, {'color'}, shimColors);
xlim( [0 tvec(end)/1e-3] );
grid(axShimwaveform, 'on');
axShimwaveform.TickLabelInterpreter='latex';
axShimwaveform.FontSize = fontSizeAx;
xticklabels( axShimwaveform, [] );
ylabel('Shim (A)', 'interpreter', 'latex');
title(['\textbf{', 'Shim', '}'], 'interpreter', 'latex');

%% Grad Pulse Diagram
[ tvec, Grad ] = makeBoxyWaveformAdjoint( wv.tvec, wv.dtvec, wv.Grad );

axGradwaveform = nexttile(4);
pGrad = plot( transpose( tvec )/1e-3,  transpose( Grad )/1e-3);
set(pGrad, 'linewidth', lw);
set(pGrad, {'color'}, gradColors);
xlim( [0 tvec(end)/1e-3] );
grid(axGradwaveform, 'on');
ylabel('Grad (mT/m)', 'interpreter', 'latex');
xlabel('Time (ms)', 'interpreter', 'latex')
title(['\textbf{', 'Grad', '}'], 'interpreter', 'latex');
legend('show', {'$G_x$','$G_y$','$G_z$'},...
    'interpreter', 'latex',...
    'location', 'northwest')
grid(axGradwaveform, 'on');
axGradwaveform.TickLabelInterpreter='latex';
axGradwaveform.FontSize = fontSizeAx;

%% Set figure window style back to default
if ~binVis
    set(groot, 'defaultfigurewindowstyle', dfws)
end


end