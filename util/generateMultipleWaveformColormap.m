function wvColors = generateMultipleWaveformColormap( numWv, rgb, ColorMax, ColorMin, OtherMax )
arguments
    numWv
    rgb
    ColorMax = 1;
    ColorMin = 0.25;
    OtherMax = 0.75;
end

if strcmpi( rgb, 'r') || strcmpi( rgb, 'red')
    rgb_idx = 1;
elseif strcmpi( rgb, 'g') || strcmpi( rgb, 'green')
    rgb_idx = 2;
elseif strcmpi( rgb, 'b') || strcmpi( rgb, 'blue')
    rgb_idx = 3;
else
    error( "Unknown color to make colormap." )
end

if mod(numWv,2) == 0
    numInt = numWv/2;
    dcDark = (ColorMax-ColorMin)/(numInt);
    interpDark = (ColorMin:dcDark:ColorMax).';
    % wvColorsDark = [...
    %     zeros(length(interpDark), 1),...
    %     zeros(length(interpDark), 1) ...
    %     interpDark,...
    %     ];
    wvColorsDark = zeros( length(interpDark), 3 );
    wvColorsDark( :, rgb_idx ) = interpDark;

    dcLight = (OtherMax-0)/(numInt-1);
    interpLight = (dcLight:dcLight:OtherMax).';
    % wvColorsLight = [...
    %     interpLight,...
    %     interpLight,...
    %     ones(length(interpLight), 1),...
    %     ];
    wvColorsLight = zeros( length(interpLight), 3 );
    wvColorsLight( :, rgb_idx ) = interpLight;

    wvColors = [ wvColorsDark; wvColorsLight ];
elseif numWv == 1
    if rgb_idx == 1
        wvColors = "#A2142F";
    elseif rgb_idx == 2
        wvColors = 	"#77AC30";
    elseif rgb_idx == 3
        wvColors = "#0072BD";
    end
else
    numInt = floor(numWv/2);
    dcDark = (ColorMax-ColorMin)/numInt;
    interpDark = (ColorMin:dcDark:ColorMax).';
    % wvColorsDark = [...
    %     zeros(length(interpDark), 1),...
    %     zeros(length(interpDark), 1),...
    %     interpDark,...
    %     ];
    wvColorsDark = zeros( length(interpDark), 3 );
    wvColorsDark( :, rgb_idx ) = interpDark;

    dcLight = (OtherMax-0)/numInt;
    interpLight = (dcLight:dcLight:OtherMax).';
    % wvColorsLight = [...
    %     interpLight,...
    %     interpLight,...
    %     ones(length(interpLight), 1),...
    %     ];
    wvColorsLight = zeros( length(interpLight), 3 );
    wvColorsLight( :, rgb_idx ) = interpLight;

    wvColors = [ wvColorsDark; wvColorsLight ];
end

wvColors = flip( wvColors, 1 );
wvColors = num2cell(wvColors, 2);

end