function [ FAFig, phaseFig ] = plotInSlice(...
    val, binVis )

%% Pre-process data

if isfield( val, "flipAngleRangePlot" )
    faRange = val.flipAngleRangePlot;
else
    faRange = [];
end

%% Determine Array Sizes
numSubj = val.numSubj;
numSlices = val.numSlices;

axSlicesFA = zeros( [ size( val.FAInPlane{ 1, 1 }, [ 1, 2] ), numSlices, numSubj ] );

if nargout > 1
    axSlicesMph = zeros( [ size( val.FAInPlane{ 1, 1 }, [1, 2] ), numSlices, numSubj ] );
end

axSlicesROI = zeros( [ size( val.FAInPlane{ 1, 1 }, [1, 2] ), numSlices, numSubj ] );

%% Get Slices
for ll = 1:numSlices
    for ss = 1:numSubj

        if val.calcInPlaneMetrics( ll, ss )

            [ ~, zCenterIdx ] = min( abs( val.zInPlane{ ll, ss } - val.sliceLocation( ll ) ) );

            axSlicesFA( :, :, ll, ss ) = squeeze( val.FAInPlane{ ll, ss }( :, :, zCenterIdx ) );
            axSlicesROI( :, :, ll, ss ) = squeeze( val.optroiInPlane{ ll, ss }( :, :, zCenterIdx ) );

            if nargout > 1
                axSlicesMph( :, :, ll, ss ) = atan2(...
                    squeeze( val.MInPlane{ ll, ss }( :, :, zCenterIdx, 2 ) ),...
                    squeeze( val.MInPlane{ ll, ss }( :, :, zCenterIdx, 1 ) ) );
            end

        else
            axSlicesFA( :, :, ll, ss ) = zeros(size( val.FAInPlane{ 1, 1 }, [ 1, 2 ] ) );
            axSlicesROI( :, :, ll, ss ) = false( size( val.FAInPlane{ 1, 1 }, [ 1, 2 ] ) );

            if nargout > 1
                axSlicesFA( :, :, ll, ss ) = zeros( size( val.FAInPlane{ 1, 1 }, [ 1, 2 ] ) );
            end
        end
    end
end

%% Plot FA
FAFig = generateFAFig( axSlicesFA, axSlicesROI, faRange, val.subjIden, binVis );

%% Plot phase
if nargout > 1
    phaseFig = generatePhFig( axSlicesMph, axSlicesROI, val.subjIden, binVis );
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function FAFig = generateFAFig( axSlicesFA, axSlicesROI, faRange, subjIden, binVis )

numSubj = size( axSlicesFA, 4 );
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

if isempty( faRange )
    maxFA = max( axSlicesFA(:),...
        [], 'all' );
    cAx = [ 0, maxFA ];
else
    cAx = [ faRange(1), faRange(2) ];
end

%% Pre plotting computation
magm = colorcet('L09');
olC = 'w';
clW = 1.75;
fsCb = 24;

figSize = [1 1 numCols*500 numRows*500];

axPercHorzReduce = 0.15;
axPercVertReduce = 0.10;

%% Initialize Figure
if binVis
    FAFig = figure('color', 'white', 'units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    FAFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

tlFA = tiledlayout( FAFig, numRows, numCols, 'tilespacing', 'compact', 'padding', 'compact' );

for ss = 1:numSubj
    
    %% Create mosaic
    axSlicesMosFASubj = rearrangeCutImgStack( squeeze( axSlicesFA( :, :, :, ss ) ), axPercHorzReduce, axPercVertReduce );
    axSlicesMosROISubj = rearrangeCutImgStack( squeeze( axSlicesROI( :, :, :, ss ) ), axPercHorzReduce, axPercVertReduce );

    %% Preprocess axes
    lastRow = ( ss > ( numCols * ( numRows - 1 ) ) );
    % firstCol = ( mod( ss-1, numCols ) == 0 );
    lastCol = ( mod( ss, numCols ) == 0 );

    if lastRow
        axFA = nexttile( tlFA, ( numCols * ( numRows - 1 ) ) + colLastRowIdx( ss - ( numCols * ( numRows - 1 ) ) ) );
    else
        axFA = nexttile( tlFA, ss );
    end

    %% Plot Magnitude Images
    hold(axFA, 'on');
    imagesc( axFA, flipud( rot90( abs(axSlicesMosFASubj), 1 ) ) );
    set( axFA, 'ydir', 'normal' );
    axis( axFA, 'image' );
    xticklabels( axFA, [] );
    yticklabels( axFA, [] );
    clim( axFA, cAx );
    colormap( axFA, magm );

    [ ~, contMagAx] = imcontour( flipud( rot90( axSlicesMosROISubj, 1 ) ) );
    contMagAx.Color = olC;
    contMagAx.LineWidth = clW;
    pause( 0.1 );

    cbMag = colorbar( axFA );
    cbMag.TickLabelInterpreter = 'latex';
    cbMag.FontSize = fsCb - 4;
    if lastCol
        ylabel( cbMag, 'Flip Angle [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
    end
    pause( 0.1 );

    title( strcat( "\textbf{", insertBefore(subjIden( ss ), "_", "\"), "}" ), 'interpreter', 'latex'  );

end

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function phaseFig = generatePhFig( axSlicesMph, axSlicesROI, subjIden, binVis )

numSubj = size( axSlicesMph, 4 );
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

%% Pre plotting computation
phasem = phasemap( 256 );
olC = 'w';
clW = 1.75;
fsCb = 24;

figSize = [1 1 numCols*500 numRows*500];

axPercHorzReduce = 0.15;
axPercVertReduce = 0.10;

%% Initialize Figure
if binVis
    phaseFig = figure('color', 'white', 'units', 'pixels',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    phaseFig = figure('color', 'white', 'Units', 'pixels',...
        'Visible','off','position', figSize);
end

tlPh = tiledlayout( phaseFig, numRows, numCols, 'tilespacing', 'compact', 'padding', 'compact' );

for ss = 1:numSubj
    
    %% Create mosaic
    axSlicesMosROISubj = rearrangeCutImgStack( squeeze( axSlicesROI( :, :, :, ss ) ), axPercHorzReduce, axPercVertReduce );
    axSlicesMosPhSubj = rearrangeCutImgStack( squeeze( axSlicesMph( :, :, :, ss ) ), axPercHorzReduce, axPercVertReduce );

    %% Preprocess axes
    lastRow = ( ss > ( numCols * ( numRows - 1 ) ) );
    % firstCol = ( mod( ss-1, numCols ) == 0 );
    lastCol = ( mod( ss, numCols ) == 0 );

    if lastRow
        axPh = nexttile( tlPh, ( numCols * ( numRows - 1 ) ) + colLastRowIdx( ss - ( numCols * ( numRows - 1 ) ) ) );
    else
        axPh = nexttile( tlPh, ss );
    end

    hold(axPh, 'on');
    imagesc( axPh, flipud( rot90( abs(axSlicesMosPhSubj), 1 ) ) );
    set( axPh, 'ydir', 'normal' );
    axis( axPh, 'image' );
    xticklabels( axPh, [] );
    yticklabels( axPh, [] );
    clim( axPh, [ -pi, pi ] );
    colormap( axPh, phasem );

    [ ~, contPhAx] = imcontour( flipud( rot90( axSlicesMosROISubj, 1 ) ) );
    contPhAx.Color = olC;
    contPhAx.LineWidth = clW;
    pause( 0.1 );
    if lastCol
        cbPh = colorbar( axPh );
        cbPh.FontSize = fsCb - 4;
        set(cbPh, 'YTick', -pi:pi/2:pi);
        set(cbPh, 'YTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        cbPh.TickLabelInterpreter = 'tex';
        ylabel( cbPh, 'Phase [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
    end
    pause( 0.1 );

    title( strcat( "\textbf{", insertBefore(subjIden( ss ), "_", "\"), "}" ), 'interpreter', 'latex'  );

end

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end
end
% ----------------------------------------------------------------------- %