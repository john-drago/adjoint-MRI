function freqFig = plotPulseFreq( opt, binVis )

%% Initialize variables
fontSizeAx = 18;

%% Get Waveforms
wv = opt.generateWaveforms( opt.pOpt, opt );

bcomp_rshp = complex( wv.breal, wv.bimag );

fvec = fvecDFT( 1/opt.dt, opt.numTimePoints, true );
bcomp_rshp_fft = fftshift( fft( bcomp_rshp, [], 2 ), 2 ) / opt.numTimePoints;

%% Define Params to Use During Plotting

numXYCoils = opt.numXYCoils;

lw = 2.5;

rfColors = num2cell(lines(numXYCoils), 2);
alphaVal = 0.6;
rfColors = cellfun(@(x) [x, alphaVal], rfColors, 'UniformOutput', false);

%%  Plot the Param Output
figSize = [ 1 1 850 400 ];

if binVis
    freqFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    freqFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

freqTL = tiledlayout( 1, 1,...
    'padding', 'compact');

%% Freq plot
axRFfreq = nexttile( freqTL, 1);
pRF = plot( transpose( fvec ) / 1e3, transpose( abs( bcomp_rshp_fft ) ) );
set(pRF, 'linewidth', lw);
set(pRF, {'color'}, rfColors);
xlim( [ min( fvec ) max( fvec ) ] / 1e3 );
grid(axRFfreq, 'on');
axRFfreq.TickLabelInterpreter='latex';
axRFfreq.FontSize = fontSizeAx;
% xticklabels( axRFfreq, [] );
xlabel('Freq (kHz)', 'interpreter', 'latex');
ylabel('DFT(RF)', 'interpreter', 'latex');
title(['\textbf{', 'DFT of RF waveforms', '}'], 'interpreter', 'latex');

%% Set figure window style back to default
if ~binVis
    set(groot, 'defaultfigurewindowstyle', dfws)
end


end