function [ db0Fig, b1pFig, brainMaskFig, bodyMaskFig ] = ...
    generateFieldMapImgs( Fields, binVis )
arguments
    Fields
    binVis = false
end

% Show dB0
rbcmap = colorcet( 'D01' );
rbcaxis = max( abs([ max(Fields.gyro/(2*pi) * Fields.DB0,[], 'all'), min(Fields.gyro/(2*pi) * Fields.DB0,[], 'all') ]) );
rbcaxis = [ -rbcaxis rbcaxis ];
pause( 0.5 );
db0Fig = plotImgStack( ( Fields.gyro/(2*pi) * Fields.DB0), rbcmap, rbcaxis, binVis );
db0CBar = findall( db0Fig, 'type', 'colorbar' );
db0CBar.Label.String = "DB0 (Hz)";

% Show B1p
lcmap = colorcet( 'L09' );
lcaxis = [ 0, max( abs(Fields.b1p)/1e-9, [], 'all' ) ];
pause( 0.5 );
b1pFig = plotImgStack( ( abs( Fields.b1p ) / 1e-9 ), lcmap, lcaxis, binVis );
b1pCBar = findall( b1pFig, 'type', 'colorbar' );
b1pCBar.Label.String = "B1+ (nT/V)";

% Show Brain Mask
binmap = gray;
binval = [0 1];
pause( 0.5 );
brainMaskFig = plotImgStack( Fields.roi_brain, binmap, binval, binVis );
brainMaskBar = findall( brainMaskFig, 'type', 'colorbar' );
brainMaskBar.Label.String = "Brain Mask";

% Shown Body Mask
pause( 0.5 );
bodyMaskFig = plotImgStack( Fields.roi_body, binmap, binval, binVis );
bodyMaskBar = findall( bodyMaskFig, 'type', 'colorbar' );
bodyMaskBar.Label.String = "Body Mask";

end


%% Helper Functions
% ----------------------------------------------------------------------- %
function imgStackFig = plotImgStack( imgStack, cmap, cval, binVis )

figSize = [ 1 1 800 600 ];

if binVis
    imgStackFig = figure('color', 'white',...
        'Visible','on');
else
    dfws = get(groot, 'defaultfigurewindowstyle');
    set(groot, 'defaultfigurewindowstyle','normal')
    imgStackFig = figure('color', 'white',...
        'Visible','off','position', figSize);
end

imgStackMosaic = rearrangeCutImgStack( imgStack, 0, 0 );

ax = axes();
imagesc( permute( imgStackMosaic, [2 1]) ); 
set(ax, 'YDir', 'normal')
axis image;
colormap(ax, cmap);
colorbar(ax);
caxis(ax, cval); %#ok

if ~binVis % change default figure window style back to original
    set(groot, 'defaultfigurewindowstyle', dfws)
end

end
% ----------------------------------------------------------------------- %
% % ----------------------------------------------------------------------- %
% function img3DFig = plot3DCutouts( imgStack, X, Y, Z, cmap, cval, binVis )
% 
% figSize = [ 1 1 800 600 ];
% 
% if binVis
%     img3DFig = figure('color', 'white',...
%         'Visible','on');
% else
%     dfws = get(groot, 'defaultfigurewindowstyle');
%     set(groot, 'defaultfigurewindowstyle','normal')
%     img3DFig = figure('color', 'white',...
%         'Visible','off','position', figSize);
% end
% 
% ax = axes();
% 
% % Rearrange the dimensions to make nice with MATLAB 3D plotting
% Xp = permute(X, [2,1,3]);
% Yp = permute(Y, [2,1,3]);
% Zp = permute(Z, [2,1,3]);
% 
% imgStackp = permute( imgStack, [2 1 3] );
% 
% hold on;
% 
% s1 = slice(Xp, Yp, Zp, imgStackp, 0, [], []);
% s1.EdgeColor = "none";
% 
% s2 = slice(Xp, Yp, Zp, imgStackp, [], 0, []);
% s2.EdgeColor = "none";
% 
% s3 = slice(Xp, Yp, Zp, imgStackp, [], [], 0);
% s3.EdgeColor = "none";
% 
% grid on
% xlabel('X Direction (m)')
% ylabel('Y Direction (m)')
% zlabel('Z Direction (m)')
% axis equal;
% 
% colormap(ax, cmap);
% % clim( db0CMap );
% colorbar(ax);
% caxis(ax, cval ); %#ok
% 
% ax.CameraPosition = [1.728817361808542   0.563928000193659   0.551258904368493];
% ax.CameraUpVector = [ 0, 0, 1 ];
% ax.CameraViewAngle = 8.3278;
% 
% if ~binVis % change default figure window style back to original
%     set(groot, 'defaultfigurewindowstyle', dfws)
% end
% 
% end
% % ----------------------------------------------------------------------- %