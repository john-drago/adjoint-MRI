function [ thruSliceFAFig, thruSlicePhFig ] = plotThruSlice( val, binVis )

%% Pre-process data
numSubj = val.numSubj;
numRows = floor( sqrt( numSubj ) );
numCols = ceil( numSubj / numRows );
numColLastRow = numSubj - ( numCols * ( numRows - 1 ) );

centerCol = numCols / 2 + 0.5;
colExtLastRow = numColLastRow / 2 ;
numDigRound = 3;
colLastRow = [...
    ceil( round( centerCol - colExtLastRow, numDigRound ) ),...
    floor( round( centerCol + colExtLastRow, numDigRound ) - 0.5 ) ];

colLastRowIdx = colLastRow( 1 ) : colLastRow( 2 );
colLastRowIdx = colLastRowIdx( colLastRowIdx ~= 0 );

%% Get colormaps
ColorMax = 1;
ColorMin = 0.25;
OtherMax = 0.75;

%% Set parameters
lw = 2;
targetlw = 4;
di_tol = 1e-8;

% Number of candidate phase windows to evaluate (1–8 supported here)
num_windows = 4;

% Font size to match database plotting
fontsize = 36;

figSize = [ 1, 1, numCols*600, numRows*300 ];

%% Initialize FA Figure
if binVis
    thruSliceFAFig = figure('color', 'white', 'units', 'pixels','Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    thruSliceFAFig = figure('color', 'white', 'Units', 'pixels','Visible','off','position', figSize);
end

phaseExtendAmt = 0.10;
boundsExtendPhase = ( [...
    val.sliceBounds( :, 1 ) - val.sliceThickness * phaseExtendAmt,...
    val.sliceBounds( :, 2 ) + val.sliceThickness * phaseExtendAmt ] ) * 1000 ;

extraBuffer = 1.0 * val.sliceThickness * 1000;
zBounds = [...
    min( boundsExtendPhase, [], 'all' ) - extraBuffer,...
    max( boundsExtendPhase, [], 'all' ) + extraBuffer ];

tl = tiledlayout( numRows, numCols, 'tilespacing', 'compact', 'padding', 'compact' );

for ss = 1:numSubj

    lastRow = ( ss > ( numCols * ( numRows - 1 ) ) );
    firstCol = ( mod( ss-1, numCols ) == 0 );

    if lastRow
        ax = nexttile( tl, ( numCols * ( numRows - 1 ) ) + colLastRowIdx( ss - ( numCols * ( numRows - 1 ) ) ) );
    else
        ax = nexttile( tl, ss );
    end

    for ll = 1:val.numSlices

        if val.calcThruMetrics( ll, ss )

            FAColors = generateMultipleWaveformColormap( size( val.xySamplePtFAVals{ ll, ss }, 2 ), 'b', ColorMax, ColorMin, OtherMax );

            FAThruSubj = val.xySamplePtFAVals{ ll, ss };
            FAtargThruSubj = val.xySamplePtFAtargVals{ ll, ss }( :, 1 );
            zThruSubj = val.xySamplePtzVals{ ll, ss }( :, 1 ) * 1000; % convert to mm

            sliceIdxMagStart = find( zThruSubj >= ( zBounds( ll, 1 ) - di_tol ), 1, 'first' );
            sliceIdxMagEnd   = find( zThruSubj <= ( zBounds( ll, 2 ) + di_tol ), 1, 'last' );

            hold( ax, 'on' ); grid( ax, 'on' );

            fap = plot( zThruSubj, FAThruSubj );
            set( fap, 'marker', 'none', 'linewidth', lw, 'LineStyle','-', {'color'}, FAColors, 'HandleVisibility','off' );
            ax.YAxis(1).Color = 'b';
            if firstCol
                ylabel( 'FA [deg]', 'interpreter', 'latex' );
            end

            px = plot( zThruSubj(sliceIdxMagStart:sliceIdxMagEnd), FAtargThruSubj(sliceIdxMagStart:sliceIdxMagEnd,1) );
            px.Color = "#D95319";
            px.Marker = "none";
            px.LineWidth = targetlw;
            px.LineStyle = ':';
            px.DisplayName = "Target";

        end
    end

    xlim( zBounds );
    title( strcat( "\textbf{", insertBefore(val.subjIden( ss ), "_", "\"), "}" ), 'interpreter', 'latex' );
    if lastRow
        xlabel( 'Thru Slice Position [mm]', 'interpreter', 'latex' );
    end

    % ----- Font sizes like database function -----
    ax.FontSize = fontsize;
    ax.XLabel.FontSize = fontsize;
    ax.YLabel.FontSize = fontsize;
    ax.Title.FontSize  = fontsize;
    if isprop(ax,'YAxis') && numel(ax.YAxis) == 2
        ax.YAxis(1).FontSize = fontsize;
        ax.YAxis(2).FontSize = fontsize;
    end
end

%% Initialize Ph Figure
if binVis
    thruSlicePhFig = figure('color', 'white', 'units', 'pixels','Visible','on');
else
    thruSlicePhFig = figure('color', 'white', 'Units', 'pixels','Visible','off','position', figSize);
end

tl = tiledlayout( numRows, numCols, 'tilespacing', 'compact', 'padding', 'compact' );

for ss = 1:numSubj

    lastRow = ( ss > ( numCols * ( numRows - 1 ) ) );
    firstCol = ( mod( ss-1, numCols ) == 0 );
    lastCol  = ( mod( ss,   numCols ) == 0 );

    if lastRow
        ax = nexttile( tl, ( numCols * ( numRows - 1 ) ) + colLastRowIdx( ss - ( numCols * ( numRows - 1 ) ) ) );
    else
        ax = nexttile( tl, ss );
    end

    for ll = 1:val.numSlices

        if val.calcThruMetrics( ll, ss )

            PhColors = generateMultipleWaveformColormap( size( val.xySamplePtFAVals{ ll, ss }, 2 ), 'r', ColorMax, ColorMin, OtherMax );

            % Compute phase and targets
            phaseMThruSubj = atan2( val.xySamplePtMyVals{ ll, ss }, val.xySamplePtMxVals{ ll, ss } );
            FAtargThruSubj = val.xySamplePtFAtargVals{ ll, ss }( :, 1 );
            zThruSubj      = val.xySamplePtzVals{ ll, ss }( :, 1 ) * 1000; % convert to mm

            % Candidate windows: each row is [lower_bound, upper_bound]
            % Note: comments avoid using the pi character; strings still use '\pi' for TeX.
            potential_windows = [
                0, 2*pi;          % [0, two*pi]
                -pi/2, 3*pi/2;    % [-pi/2, 3pi/2]
                -pi, pi;          % [-pi, pi]
                -3*pi/2, pi/2     % [-3pi/2, pi/2]
            ];
            windows_to_evaluate = potential_windows(1:min(num_windows, size(potential_windows,1)), :);

            % Indices for plotting bounds
            sliceIdxPhStart = find( zThruSubj >= ( boundsExtendPhase( ll, 1 ) - di_tol ), 1, 'first' );
            sliceIdxPhEnd   = find( zThruSubj <= ( boundsExtendPhase( ll, 2 ) + di_tol ), 1, 'last' );
            sliceIdxMagStart = find( zThruSubj >= ( zBounds( ll, 1 ) - di_tol ), 1, 'first' );
            sliceIdxMagEnd   = find( zThruSubj <= ( zBounds( ll, 2 ) + di_tol ), 1, 'last' );

            % For jump cost, only evaluate inside actual slice bounds (mm)
            sliceIdxJumpStart = find( zThruSubj >= ( val.sliceBounds(ll,1)*1000 - di_tol ), 1, 'first' );
            sliceIdxJumpEnd   = find( zThruSubj <= ( val.sliceBounds(ll,2)*1000 + di_tol ), 1, 'last' );
            phases_slice_for_jumps = phaseMThruSubj(sliceIdxJumpStart:sliceIdxJumpEnd, :);

            % Choose window that minimizes total phase jumps along z
            best_window_idx = 1; min_jump_cost = inf;
            for w = 1:size(windows_to_evaluate,1)
                window_bounds = windows_to_evaluate(w,:);
                phases_in_window = convertToWindow(phases_slice_for_jumps, window_bounds);
                jump_cost = calculateJumpCost(phases_in_window);
                if jump_cost < min_jump_cost
                    min_jump_cost = jump_cost;
                    best_window_idx = w;
                end
            end
            chosen_window = windows_to_evaluate(best_window_idx,:);

            % Convert full phase array to chosen window
            phaseMThruSubj = convertToWindow(phaseMThruSubj, chosen_window);

            % Tick setup based on chosen window
            ylim_bounds = [chosen_window(1) - 0.05*pi, chosen_window(2) + 0.05*pi];
            tick_spacing = pi/2;
            yticks = ceil(chosen_window(1)/tick_spacing)*tick_spacing : tick_spacing : ...
                     floor(chosen_window(2)/tick_spacing)*tick_spacing;
            yticklabels = arrayfun(@(v) formatPiLabel(v), yticks, 'UniformOutput', false);

            % Phase (left axis)
            yyaxis(ax,'left'); hold(ax,'on'); grid(ax,'on');
            ylim(ylim_bounds);
            set(ax,'YTick',yticks,'YTickLabel',yticklabels);
            pp = plot( zThruSubj(sliceIdxPhStart:sliceIdxPhEnd), ...
                       phaseMThruSubj(sliceIdxPhStart:sliceIdxPhEnd,:) );
            set(pp, 'marker','none', {'color'}, PhColors, 'LineWidth', lw, 'LineStyle','-', 'HandleVisibility','off');
            ax.YAxis(1).Color = '#A2142F';
            ax.TickLabelInterpreter = 'tex';
            if firstCol
                ylabel( ax, 'In-plane phase [rad]', 'interpreter', 'latex' );
            end

            % FA target (right axis)
            yyaxis(ax,'right'); hold(ax,'on'); grid(ax,'on');
            px = plot( zThruSubj(sliceIdxMagStart:sliceIdxMagEnd), ...
                       FAtargThruSubj(sliceIdxMagStart:sliceIdxMagEnd,1) );
            px.Color = "b"; px.Marker = "none"; px.LineWidth = targetlw; px.LineStyle = ':';
            ylim( 1.05 * [ 0, max( FAtargThruSubj(sliceIdxMagStart:sliceIdxMagEnd,1), [], "all" ) ] );
            ax.YAxis(2).Color = 'b';
            if lastCol
                ylabel( 'FA [deg]', 'interpreter', 'latex' );
            end

        end
    end

    title( strcat( "\textbf{", insertBefore(val.subjIden( ss ), "_", "\"), "}" ), 'interpreter', 'latex' );
    xlim( zBounds );
    if lastRow
        xlabel( 'Thru Slice Position [mm]', 'interpreter', 'latex' );
    end

    ax.FontSize = fontsize;
    ax.XLabel.FontSize = fontsize;
    ax.YLabel.FontSize = fontsize;
    ax.Title.FontSize  = fontsize;
    if isprop(ax,'YAxis') && numel(ax.YAxis) == 2
        ax.YAxis(1).FontSize = fontsize;
        ax.YAxis(2).FontSize = fontsize;
    end
end

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end

%% Helper Functions

% ----------------------------------------------------------------------- %
function phases_out = convertToWindow(phases_in, window_bounds)
    % Convert phases to the specified window [lower_bound, upper_bound]
    lower_bound = window_bounds(1);
    upper_bound = window_bounds(2);
    window_width = upper_bound - lower_bound;

    % Shift phases so lower_bound maps to zero
    phases_shifted = phases_in - lower_bound;

    % Wrap to [0, window_width]
    phases_wrapped = mod(phases_shifted, window_width);

    % Shift back to [lower_bound, upper_bound]
    phases_out = phases_wrapped + lower_bound;
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function cost = calculateJumpCost(phases)
    % Total absolute jump between consecutive z-samples (summed over columns)
    if size(phases,1) < 2
        cost = 0;
        return;
    end
    diffs = diff(phases,1,1);
    cost = sum(abs(diffs), 'all');
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function label = formatPiLabel(value)
    % Create TeX label as multiples of \pi
    if abs(value) < 1e-10
        label = '0';
        return;
    end
    pi_multiple = value / pi;

    if abs(pi_multiple - round(pi_multiple)) < 1e-10
        k = round(pi_multiple);
        if k == 1
            label = '\pi';
        elseif k == -1
            label = '-\pi';
        else
            label = sprintf('%d\\pi', k);
        end
        return;
    end

    common = [0.5, 1.5, -0.5, -1.5];
    labels = {'\pi/2','3\pi/2','-\pi/2','-3\pi/2'};
    for ii = 1:numel(common)
        if abs(pi_multiple - common(ii)) < 1e-10
            label = labels{ii};
            return;
        end
    end

    [num, den] = rat(pi_multiple, 1e-3);
    if den == 1
        if num == 1,  label = '\pi';
        elseif num == -1, label = '-\pi';
        else, label = sprintf('%d\\pi', num);
        end
    else
        if num == 1,      label = sprintf('\\pi/%d', den);
        elseif num == -1, label = sprintf('-\\pi/%d', den);
        else,             label = sprintf('%d\\pi/%d', num, den);
        end
    end
end
% ----------------------------------------------------------------------- %
