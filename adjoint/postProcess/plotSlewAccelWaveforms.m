function slewAccelFig = plotSlewAccelWaveforms( opt, binVis )

%% Initialize variables
fontSizeAx = 12;

%% Get Waveforms
wv = opt.generateWaveforms( opt.pOpt, opt );
wv = calculateWaveformSlewAccel( wv );

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
figSize = [ 1 1 2000 800 ];

if binVis
    slewAccelFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    slewAccelFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

PSDTL = tiledlayout( 4, 2,...
    'padding', 'compact');

%% brealphasor slew plot
[ tvec, brealphasor_slew ] = makeBoxyWaveformAdjoint( wv.tvec_slew, wv.dtvec_slew, wv.brealphasor_slew );

axRFbreal = nexttile(PSDTL, 1);
pRF = plot( transpose( tvec )/1e-3, transpose( brealphasor_slew ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsReal);
xlim( [0 tvec(end)/1e-3] );
ylim( [-max(abs(brealphasor_slew), [], 'all'), max(abs(brealphasor_slew), [], 'all')] );
grid(axRFbreal, 'on');
axRFbreal.TickLabelInterpreter='latex';
axRFbreal.FontSize = fontSizeAx;
xticklabels( axRFbreal, [] );
ylabel('Slew $b_{\rm real}$ (V/sec)', 'interpreter', 'latex');
title(['\textbf{', 'Slew RF real', '}'], 'interpreter', 'latex');

%% brealphasor acc plot
[ tvec, brealphasor_accel ] = makeBoxyWaveformAdjoint( wv.tvec_accel, wv.dtvec_accel, wv.brealphasor_accel );

axRFbreal = nexttile(PSDTL, 2);
pRF = plot( transpose( tvec )/1e-3, transpose( brealphasor_accel ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsReal);
xlim( [0 tvec(end)/1e-3] );
ylim( [-max(abs(brealphasor_accel), [], 'all'), max(abs(brealphasor_accel), [], 'all')] );
grid(axRFbreal, 'on');
axRFbreal.TickLabelInterpreter='latex';
axRFbreal.FontSize = fontSizeAx;
xticklabels( axRFbreal, [] );
ylabel('Accel $b_{\rm real}$ (V/sec$^2$)', 'interpreter', 'latex');
title(['\textbf{', 'Accel RF real', '}'], 'interpreter', 'latex');

%% bimagphasor slew plot
[ tvec, bimagphasor_slew ] = makeBoxyWaveformAdjoint( wv.tvec_slew, wv.dtvec_slew, wv.bimagphasor_slew );

axRFbimag = nexttile(PSDTL, 3);
pRF = plot( transpose( tvec )/1e-3, transpose( bimagphasor_slew ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsImag);
xlim( [0 tvec(end)/1e-3] );
ylim( [-max(abs(bimagphasor_slew), [], 'all'), max(abs(bimagphasor_slew), [], 'all')] );
grid(axRFbimag, 'on');
axRFbimag.TickLabelInterpreter='latex';
axRFbimag.FontSize = fontSizeAx;
xticklabels( axRFbimag, [] );
ylabel('Slew $b_{\rm imag}$ (V/sec)', 'interpreter', 'latex');
title(['\textbf{', 'Slew RF imag', '}'], 'interpreter', 'latex');

%% bimagphasor acc plot
[ tvec, bimagphasor_accel ] = makeBoxyWaveformAdjoint( wv.tvec_accel, wv.dtvec_accel, wv.bimagphasor_accel );

axRFbimag = nexttile(PSDTL, 4);
pRF = plot( transpose( tvec )/1e-3, transpose( bimagphasor_accel ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColorsImag);
xlim( [0 tvec(end)/1e-3] );
ylim( [-max(abs(bimagphasor_accel), [], 'all'), max(abs(bimagphasor_accel), [], 'all')] );
grid(axRFbimag, 'on');
axRFbimag.TickLabelInterpreter='latex';
axRFbimag.FontSize = fontSizeAx;
xticklabels( axRFbimag, [] );
ylabel('Accel $b_{\rm imag}$ (V/sec$^2$)', 'interpreter', 'latex');
title(['\textbf{', 'Accel RF imag', '}'], 'interpreter', 'latex');

%% Shim Slew plot
[ tvec, Shim_slew ] = makeBoxyWaveformAdjoint( wv.tvec_slew, wv.dtvec_slew, wv.Shim_slew );

axShimwaveform = nexttile(PSDTL, 5);
pShim = plot( transpose( tvec )/1e-3,  transpose( Shim_slew ) );
set(pShim, 'linewidth', lw);
set(pShim, {'color'}, shimColors);
xlim( [0 tvec(end)/1e-3] );
% if all( abs(Shim_slew) == 0, 'all' )  
%     ylim( [ -1, 1 ] )
% else
%     ylim( [-max(abs(Shim_slew), [], 'all'), max(abs(Shim_slew), [], 'all')] );
% end
grid(axShimwaveform, 'on');
axShimwaveform.TickLabelInterpreter='latex';
axShimwaveform.FontSize = fontSizeAx;
xticklabels( axShimwaveform, [] );
ylabel('Slew Shim (A/sec)', 'interpreter', 'latex');
title(['\textbf{', 'Slew Shim', '}'], 'interpreter', 'latex');

%% Shim Acc plot
[ tvec, Shim_accel ] = makeBoxyWaveformAdjoint( wv.tvec_accel, wv.dtvec_accel, wv.Shim_accel );

axShimwaveform = nexttile(PSDTL, 6);
pShim = plot( transpose( tvec )/1e-3,  transpose( Shim_accel ));
set(pShim, 'linewidth', lw);
set(pShim, {'color'}, shimColors);
xlim( [0 tvec(end)/1e-3] );
% ylim( [-max(abs(Shim_accel), [], 'all'), max(abs(Shim_accel), [], 'all')] );
grid(axShimwaveform, 'on');
axShimwaveform.TickLabelInterpreter='latex';
axShimwaveform.FontSize = fontSizeAx;
xticklabels( axShimwaveform, [] );
ylabel('Accel Shim (A/sec$^2$)', 'interpreter', 'latex');
title(['\textbf{', 'Accel Shim', '}'], 'interpreter', 'latex');

%% Grad Slew Plot
[ tvec, Grad_slew ] = makeBoxyWaveformAdjoint( wv.tvec_slew, wv.dtvec_slew, wv.Grad_slew );

axGradwaveform = nexttile(PSDTL, 7);
pGrad = plot( transpose( tvec )/1e-3,  transpose( Grad_slew ) );
set(pGrad, 'linewidth', lw);
set(pGrad, {'color'}, gradColors);
xlim( [0 tvec(end)/1e-3] );
% ylim( [-max(abs(Grad_slew), [], 'all'), max(abs(Grad_slew), [], 'all')] );
grid(axGradwaveform, 'on');
ylabel('Slew Grad (T/(m $\cdot$ sec))', 'interpreter', 'latex');
xlabel('Time (ms)', 'interpreter', 'latex')
title(['\textbf{', 'Slew Grad', '}'], 'interpreter', 'latex');
legend('show', {'$G_x$','$G_y$','$G_z$'},...
    'interpreter', 'latex',...
    'location', 'northwest')
grid(axGradwaveform, 'on');
axGradwaveform.TickLabelInterpreter='latex';
axGradwaveform.FontSize = fontSizeAx;

%% Grad Acc Plot
[ tvec, Grad_accel ] = makeBoxyWaveformAdjoint( wv.tvec_accel, wv.dtvec_accel, wv.Grad_accel );

axGradwaveform = nexttile(PSDTL, 8);
pGrad = plot( transpose( tvec )/1e-3,  transpose( Grad_accel ) );
set(pGrad, 'linewidth', lw);
set(pGrad, {'color'}, gradColors);
xlim( [0 tvec(end)/1e-3] );
% ylim( [-max(abs(Grad_accel), [], 'all'), max(abs(Grad_accel), [], 'all')] );
grid(axGradwaveform, 'on');
ylabel('Accel Grad (T/(m $\cdot$ sec$^2$))', 'interpreter', 'latex');
xlabel('Time (ms)', 'interpreter', 'latex')
title(['\textbf{', 'Accel Grad', '}'], 'interpreter', 'latex');
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