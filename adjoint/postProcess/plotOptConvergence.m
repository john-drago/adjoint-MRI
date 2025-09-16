function convFig = plotOptConvergence( opt, binVis )

%% Initialize plotting
prevColor = "#77AC30";
prevShape = "d";

lw = 1.5;
faceAlpha = 0.5;

markerSize = 50;

ylimMinInitGuess = 0.025;
ylimMinInitData = 0.975 * min( opt.fval_conv_iters( : ) );
if ylimMinInitData < ylimMinInitGuess
    ylimMin = max( [ ylimMinInitData, 0 ] );
else
    ylimMin = ylimMinInitData;
end

ylimMaxInitGuess = 0.4;
ylimMaxInitData = 1.025 * max( opt.fval_conv_iters( : ) );

ylimMax = max(ylimMaxInitGuess, ylimMaxInitData);

% if ylimMaxInitData > ylimMaxInitGuess
%     ylimMax = min( [ ylimMaxInitData, 1 ] );
% else
%     ylimMax = ylimMaxInitData;
% end

%% Initialize figure
figSize = [ 1 1 1500 550 ];
if binVis
    convFig = figure('color', 'white', 'units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    convFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

tiledlayout( 1, 2, "TileSpacing", "compact", "Padding", "compact" );

%% fval Versus Wall Time
tiax = nexttile;
hold( tiax, "on" );
tiax.XScale  = "linear";
tiax.YScale  = "linear";

prevsc = scatter( opt.optTime_conv_iters, opt.fval_conv_iters, markerSize );
prevsc.Marker = prevShape;
prevsc.LineWidth = lw;
prevsc.MarkerEdgeColor = prevColor;
prevsc.MarkerFaceColor = prevColor;
prevsc.MarkerFaceAlpha = faceAlpha;

ylabel( "Objective Function", "Interpreter", "latex" );
xlabel( "Optimization Time (secs)", "Interpreter", "latex" );

ylim( [ ylimMin, ylimMax ] );

grid( tiax, "on" );

%% fval Versus Iteration
tiax = nexttile;
hold( tiax, "on" );
tiax.XScale  = "linear";
tiax.YScale  = "linear";

prevsc = scatter( opt.funccount_conv_iters, opt.fval_conv_iters, markerSize );
prevsc.Marker = prevShape;
prevsc.LineWidth = lw;
prevsc.MarkerEdgeColor = prevColor;
prevsc.MarkerFaceColor = prevColor;
prevsc.MarkerFaceAlpha = faceAlpha;

xlabel( "Function Count", "Interpreter", "latex" );
ylim( [ ylimMin, ylimMax ] );

grid( tiax, "on" );


%% Put title
% fontSize = 26;
% title( tl, "\textbf{Optimization Convergence Plot}",...
%     "fontsize", fontSize,  "interpreter", "latex" );

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %