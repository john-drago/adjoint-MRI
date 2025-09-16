function [ wvforminterp, tvecinterp ] = interpolateWaveforms( tvec, wvform, tvecinterp, method, extrap )
% This function will interpolate waveforms at the times specified by
% tvecinterp using the wvform array and the corresponding tvec pattern.

if nargin < 4
    method = 'pchip';
end
if nargin < 5
    extrap = nan;
end

tvec = tvec( : );
tvecinterp = tvecinterp( : );

if size( wvform, 1 ) == length( tvec )
elseif size( wvform, 2 ) == length( tvec )
    wvform = transpose( wvform );
else
    error( "Size of waveform not compatible for interpolation." )
end

wvforminterp = ( interp1( tvec, wvform, tvecinterp, method, extrap ) ).';

% % determine number of waveforms
% numWv = size( wvform, 1 );
% 
% wvforminterp = zeros( numWv, length( tvecinterp ) );
% 
% for nn = 1:numWv
%     wvforminterp( nn, : ) = interp1( tvec, wvform( nn, : ), tvecinterp );
% end

end