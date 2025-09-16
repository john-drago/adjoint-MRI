function psdFig = plotAdjointPSD( val, binVis )

%% Initialize variables
fontSizeAx = 12;

%% Get Waveforms
tvec = val.wv.tvec;
Grad = val.wv.Grad;
Shim = val.wv.Shim;
RFphasor = val.wv.RFphasor;
Freq = val.wv.Freq;
numZCoils = size( Shim , 1);
numXYCoils = size( RFphasor, 1 );
tend = tvec(end);

%% Define Params to Use During Plotting
lw = 2.5;

ColorMax = 1;
ColorMin = 0.25;
OtherMax = 0.75;

shimColors = generateMultipleWaveformColormap( numZCoils, 'b', ColorMax, ColorMin, OtherMax );
rfMagColors = generateMultipleWaveformColormap( numXYCoils, 'g', ColorMax, ColorMin, OtherMax );
rfPhColors = generateMultipleWaveformColormap( numXYCoils, 'r', ColorMax, ColorMin, OtherMax );

gradColors = num2cell(lines(4), 2);
gradColors = gradColors(2:end);

%%  Plot the Param Output
figSize = [ 1 1 1100 625 ];

if binVis
    psdFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    psdFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

PSDTL = tiledlayout( 5, 1,...
    'padding', 'compact');

%% RF Pulse Diagram
axRFwaveform = nexttile(PSDTL, 1);
pRF = plot( tvec/1e-3,  abs(RFphasor));
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfMagColors);
xlim( [0 tend/1e-3] );
ylim( [0 max(abs(abs(RFphasor)), [], 'all')] );
grid(axRFwaveform, 'on');
axRFwaveform.TickLabelInterpreter='latex';
axRFwaveform.FontSize = fontSizeAx;
xticklabels( axRFwaveform, [] );
ylabel('$|{\rm RF}|$ (V)', 'interpreter', 'latex');
title(['\textbf{', 'Pulse Sequence Diagram', '}'], 'interpreter', 'latex');

%% RF Pulse Diagram
axRFwaveform = nexttile(PSDTL, 2);
pRF = plot( tvec/1e-3,  angle(RFphasor));
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfPhColors);
xlim( [0 tend/1e-3] );
ylim( [-pi, pi] * 1.10 );
grid(axRFwaveform, 'on');
axRFwaveform.TickLabelInterpreter='latex';
axRFwaveform.FontSize = fontSizeAx;
xticklabels( axRFwaveform, [] );
ylabel('$\angle {\rm RF}$ (rad)', 'interpreter', 'latex');

%% RF Freq Diagram
axFreqwaveform = nexttile(PSDTL, 3);
pRF = plot( tvec/1e-3,  Freq/1e3);
set(pRF, 'linewidth', lw);
set(pRF, 'color', '#A2142F');
xlim( [0 tend/1e-3] );
grid(axFreqwaveform, 'on');
axFreqwaveform.TickLabelInterpreter='latex';
axFreqwaveform.FontSize = fontSizeAx;
xticklabels( axFreqwaveform, [] );
ylabel('$\Delta f_{xy}$ (kHz)', 'interpreter', 'latex');

%% Shim Pulse Diagram
axShimwaveform = nexttile(PSDTL, 4);
pShim = plot( tvec/1e-3,  Shim);
set(pShim, 'linewidth', lw);
set(pShim, {'color'}, shimColors);
xlim( [0 tend/1e-3] );
grid(axShimwaveform, 'on');
axShimwaveform.TickLabelInterpreter='latex';
axShimwaveform.FontSize = fontSizeAx;
xticklabels( axShimwaveform, [] );
ylabel('Shim (A)', 'interpreter', 'latex');

%% Grad Pulse Diagram
axGradwaveform = nexttile( PSDTL, 5);
pGrad = plot( tvec/1e-3,  Grad/1e-3);
set(pGrad, 'linewidth', lw);
set(pGrad, {'color'}, gradColors);
xlim( [0 tend/1e-3] );
grid(axGradwaveform, 'on');
ylabel('Grad (mT/m)', 'interpreter', 'latex');
xlabel('Time (ms)', 'interpreter', 'latex')
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