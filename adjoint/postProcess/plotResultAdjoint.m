function [ FAFig, phaseFig ] = plotResultAdjoint(...
    val, pulse, pulsetype, binVis )

%% Get database information
opt_roi_db = val.opt_roi_array;
zi_db = val.zi;

if strcmpi( pulsetype, "opt" )
    FA_db = val.FA_opt_array;
elseif strcmpi( pulsetype, "BCHP" )
    FA_db = val.FA_BCHP_array;
else
    error( "Unknown 'pulsetype'." )
end

if nargout > 1
    if strcmpi( pulsetype, "opt" )
        Mopt_db = val.M_opt_array;
    elseif strcmpi( pulsetype, "BCHP" )
        Mopt_db = val.M_BCHP_array;
    else
        error( "Unknown 'pulsetype'." )
    end
end

if isfield( val, "flipAngleRangePlot" )
    faRange = val.flipAngleRangePlot;
else
    faRange = [];
end

%% Determine Array Sizes
numSubj = size( opt_roi_db, 4 );
ziAll = cell2mat( zi_db(:) );

ziMin = max( ( min(ziAll) - 2 ), 1 );
ziMax = min( ( max(ziAll) + 2 ), size(FA_db, 3) );
zi = ziMin : 1 : ziMax;

axSlicesFA = zeros( [ size( FA_db, [1 2] ), numSubj ] );
coSlicesFA = zeros( [ size( FA_db, 1 ), length(zi), numSubj ] );
saSlicesFA = zeros( [ size( FA_db, 2 ), length(zi), numSubj ] );

if nargout > 1
    axSlicesMph = zeros( [ size( Mopt_db, [1 2] ), numSubj ] );
    coSlicesMph = zeros( [ size( Mopt_db, 1 ), length(zi), numSubj ] );
    saSlicesMph = zeros( [ size( Mopt_db, 2 ), length(zi), numSubj ] );
end

axSlicesROI = zeros( [ size( opt_roi_db, [1 2] ), numSubj ] );
coSlicesROI = zeros( [ size( opt_roi_db, 1 ), length(zi), numSubj ] );
saSlicesROI = zeros( [ size( opt_roi_db, 2 ), length(zi), numSubj ] );

%% Get Slices
for nn = 1:numSubj
    xcenter = ceil( length(val.x)/2 );
    ycenter = ceil( length(val.y)/2 );

    if pulse.sliceSelective
        sliceLocsIdx = zeros( pulse.numSlices, 1 );
        for ss = 1:pulse.numSlices
            [ ~, sliceLocsIdx( ss ) ] = min( abs( val.z - pulse.sliceLocation( ss ) ) );
        end
        zcenter = sliceLocsIdx( floor( mean( 1:length( sliceLocsIdx ) ) ) );
    else
        zcenter = ceil( median( find( squeeze( sum( opt_roi_db( :, :, :, nn ), [ 1 2 ] ) ) ) ) );
    end
    
    axSliceFA = squeeze( FA_db( :, :, zcenter, nn ) );
    axSlicesFA( :, :, nn ) = axSliceFA;
    axSliceROI = squeeze( opt_roi_db( :, :, zcenter, nn ) );
    axSlicesROI( :, :, nn ) = axSliceROI;

    coSliceFA = squeeze( FA_db( :, ycenter, :, nn ) );
    coSlicesFA( :, :, nn ) = coSliceFA( :, zi );
    coSliceROI = squeeze( opt_roi_db( :, ycenter, :, nn ) );
    coSlicesROI( :, :, nn ) = coSliceROI( :, zi );

    saSliceFA = squeeze( FA_db( xcenter, :, :, nn ) );
    saSlicesFA( :, :, nn ) = saSliceFA( :, zi );
    saSliceROI = squeeze( opt_roi_db( xcenter, :, :, nn ) );
    saSlicesROI( :, :, nn ) = saSliceROI( :, zi );
    
    if nargout > 1
        axSliceMph = atan2(...
            squeeze( Mopt_db( :, :, zcenter, 2, nn ) ),...
            squeeze( Mopt_db( :, :, zcenter, 1, nn ) ) );
        coSliceMph = atan2(...
            squeeze( Mopt_db( :, ycenter, :, 2, nn ) ),...
            squeeze( Mopt_db( :, ycenter, :, 1, nn ) ) );
        saSliceMph = atan2(...
            squeeze( Mopt_db( xcenter, :, :, 2, nn ) ),...
            squeeze( Mopt_db( xcenter, :, :, 1, nn ) ) );

        axSlicesMph( :, :, nn ) = axSliceMph;
        coSlicesMph( :, :, nn ) = coSliceMph( :, zi );
        saSlicesMph( :, :, nn ) = saSliceMph( :, zi );
    end

end

%% Create Mosaic
axPercHorzReduce = 0.15;
axPercVertReduce = 0.10;

coPercHorzReduce = 0.15;
coPercVertReduce = 0.00;

saPercHorzReduce = 0.07;
saPercVertReduce = 0.00;

axSlicesMosFA = rearrangeCutImgStack( axSlicesFA, axPercHorzReduce, axPercVertReduce );
coSlicesMosFA = rearrangeCutImgStack( coSlicesFA, coPercHorzReduce, coPercVertReduce );
saSlicesMosFA = rearrangeCutImgStack( saSlicesFA, saPercHorzReduce, saPercVertReduce );

axSlicesROIMos = rearrangeCutImgStack( axSlicesROI, axPercHorzReduce, axPercVertReduce );
coSlicesROIMos = rearrangeCutImgStack( coSlicesROI, coPercHorzReduce, coPercVertReduce );
saSlicesROIMos = rearrangeCutImgStack( saSlicesROI, saPercHorzReduce, saPercVertReduce );

if nargout > 1
    axSlicesMosMph = rearrangeCutImgStack( axSlicesMph, axPercHorzReduce, axPercVertReduce );
    coSlicesMosMph = rearrangeCutImgStack( coSlicesMph, coPercHorzReduce, coPercVertReduce );
    saSlicesMosMph = rearrangeCutImgStack( saSlicesMph, saPercHorzReduce, saPercVertReduce );
end

%% Generate Flip Angle Fig
FAFig = generateFAFig( axSlicesMosFA, coSlicesMosFA, saSlicesMosFA,...
    axSlicesROIMos, coSlicesROIMos, saSlicesROIMos, faRange, binVis );

%% Generate Phase Fig
if nargout > 1
    phaseFig = generatePhaseFig( axSlicesMosMph, coSlicesMosMph, saSlicesMosMph,...
        axSlicesROIMos, coSlicesROIMos, saSlicesROIMos,...
        binVis );
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function phaseFig = generatePhaseFig( axSlicesMos, coSlicesMos, saSlicesMos,...
    axSlicesROIMos, coSlicesROIMos, saSlicesROIMos,...
    binVis )

%% Determine phase bounds
cAx = [ -pi, pi ];

%% Pre plotting computation
phasem = phasemap( 256 );
olC = 'w';
clW = 1.75;
fsCb = 24;

figSize = [1 1 500 800];

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

tl = tiledlayout( 3, 1, 'tilespacing', 'compact', 'padding', 'compact' );

%% Plot Axial Images
ax1 = nexttile( tl, 1 );
hold(ax1, 'on');
imagesc( ax1, flipud( rot90( abs(axSlicesMos), 1 ) ) );
set( ax1, 'ydir', 'normal' );
axis( ax1, 'image' );
xticklabels( ax1, [] );
yticklabels( ax1, [] );
clim( ax1, cAx );
colormap( ax1, phasem );

[ ~, contMagAx] = imcontour( flipud( rot90( axSlicesROIMos, 1 ) ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbPh = colorbar( ax1 );
cbPh.TickLabelInterpreter = 'latex';
cbPh.FontSize = fsCb - 4;
ylabel( cbPh, 'Phase [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
set(cbPh, 'YTick', -pi:pi/2:pi);
set(cbPh, 'YTickLabel', {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'});
pause( 0.1 );

%% Plot Sagittal Images
ax2 = nexttile( tl, 2 );
hold(ax2, 'on');
imagesc( ax2, fliplr( rot90( abs(saSlicesMos), 1 ) ) );
set( ax2, 'ydir', 'normal' );
axis( ax2, 'image' );
xticklabels( ax2, [] );
yticklabels( ax2, [] );
clim( ax2, cAx );
colormap( ax2, phasem );

[ ~, contMagAx] = imcontour( fliplr( rot90( saSlicesROIMos, 1 ) ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbPh = colorbar( ax2 );
cbPh.TickLabelInterpreter = 'latex';
cbPh.FontSize = fsCb - 4;
ylabel( cbPh, 'Phase [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
set(cbPh, 'YTick', -pi:pi/2:pi);
set(cbPh, 'YTickLabel', {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'});
pause( 0.1 );

%% Plot Coronal Images
ax3 = nexttile( tl, 3 );
hold(ax3, 'on');
imagesc( ax3, rot90( abs(coSlicesMos), 1 ) );
set( ax3, 'ydir', 'normal' );
axis( ax3, 'image' );
xticklabels( ax3, [] );
yticklabels( ax3, [] );
clim( ax3, cAx );
colormap( ax3, phasem );

[ ~, contMagAx] = imcontour( rot90( coSlicesROIMos, 1 ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbPh = colorbar( ax3 );
cbPh.TickLabelInterpreter = 'latex';
cbPh.FontSize = fsCb - 4;
ylabel( cbPh, 'Phase [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
set(cbPh, 'YTick', -pi:pi/2:pi);
set(cbPh, 'YTickLabel', {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'});
pause( 0.1 );

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function FAFig = generateFAFig( axSlicesMos, coSlicesMos, saSlicesMos,...
    axSlicesROIMos, coSlicesROIMos, saSlicesROIMos,...
    faRange, binVis )

%% Determine flip angle bounds
if isempty(faRange)
    maxFA = max( [ axSlicesMos(:); coSlicesMos(:); saSlicesMos(:) ],...
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

figSize = [1 1 500 800];

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

tl = tiledlayout( 3, 1, 'tilespacing', 'compact', 'padding', 'compact' );

%% Plot Axial Images
ax1 = nexttile( tl, 1 );
hold(ax1, 'on');
imagesc( ax1, flipud( rot90( abs(axSlicesMos), 1 ) ) );
set( ax1, 'ydir', 'normal' );
axis( ax1, 'image' );
xticklabels( ax1, [] );
yticklabels( ax1, [] );
clim( ax1, cAx );
colormap( ax1, magm );

[ ~, contMagAx] = imcontour( flipud( rot90( axSlicesROIMos, 1 ) ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbMag = colorbar( ax1 );
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = fsCb - 4;
ylabel( cbMag, 'Flip Angle [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
pause( 0.1 );

%% Plot Sagittal Images
ax2 = nexttile( tl, 2 );
hold(ax2, 'on');
imagesc( ax2, fliplr( rot90( abs(saSlicesMos), 1 ) ) );
set( ax2, 'ydir', 'normal' );
axis( ax2, 'image' );
xticklabels( ax2, [] );
yticklabels( ax2, [] );
clim( ax2, cAx );
colormap( ax2, magm );

[ ~, contMagAx] = imcontour( fliplr( rot90( saSlicesROIMos, 1 ) ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbMag = colorbar( ax2 );
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = fsCb - 4;
ylabel( cbMag, 'Flip Angle [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
pause( 0.1 );

%% Plot Coronal Images
ax3 = nexttile( tl, 3 );
hold(ax3, 'on');
imagesc( ax3, rot90( abs(coSlicesMos), 1 ) );
set( ax3, 'ydir', 'normal' );
axis( ax3, 'image' );
xticklabels( ax3, [] );
yticklabels( ax3, [] );
clim( ax3, cAx );
colormap( ax3, magm );

[ ~, contMagAx] = imcontour( rot90( coSlicesROIMos, 1 ) );
contMagAx.Color = olC;
contMagAx.LineWidth = clW;
pause( 0.1 );

cbMag = colorbar( ax3 );
cbMag.TickLabelInterpreter = 'latex';
cbMag.FontSize = fsCb - 4;
ylabel( cbMag, 'Flip Angle [deg]', 'interpreter', 'latex', 'fontsize', fsCb );
pause( 0.1 );

%% End Processes
if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %