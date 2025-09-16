function [ wvform_acc, dtvec_acc, tvec_acc ] = calculateAccel( dtvec, tvec, wvform )
% This function will calculate the acceleration of a set of waveforms, where
% each row is a distinct channel's waveform.


% calculate forward difference
dim = 2;
n = 1;


%% Reshape arrays
dtvec = transpose( dtvec( : ) );
tvec = transpose( tvec( : ) );
if size( wvform, 2 ) == length( tvec )
elseif size( wvform, 1 ) == length( tvec )
    wvform = transpose( wvform );
else
    error( "Size of waveform not compatible for interpolation." )
end

%% Slew Calculation
% expand wvform, tvec and dtvec
wvform_slew_expand = [ zeros( size( wvform, 1 ), 1 ), wvform, zeros( size( wvform, 1 ), 1 ) ];
tvec_slew_expand = [ 0, tvec, tvec(end)+dtvec(end)/2 ];

dtvec_slew = diff( tvec_slew_expand, n, dim );

tvec_slew = [ 0, tvec(1:(end))] + dtvec_slew / 2;
% tvec_slew = [ 0, cumsum( dtvec_slew(1:(end-1)) )] + dtvec_slew / 2;

tvec_slew_div = tvec_slew;
tvec_slew_div( 1 ) = tvec_slew_expand( 1 );
tvec_slew_div( end ) = tvec_slew_expand( end );

dtvec_slew_div = dtvec_slew;
dtvec_slew_div( 1 ) = dtvec( 1 );
dtvec_slew_div( end ) = dtvec( end );
wvform_slew = diff( wvform_slew_expand, n, dim ) ./ dtvec_slew_div;

% assign zeros where we divide by 0
wvform_slew( isnan( wvform_slew ) ) = 0;

%% Acceleration Calculation
dtvec_acc = diff( tvec_slew_div, n, dim );
wvform_acc = diff( wvform_slew, n, dim ) ./ dtvec_acc;
tvec_acc = tvec_slew_div(1:(end-1)) + dtvec_acc/2;

% [ wvform_slew, dtvec_slew, tvec_slew ] = calculateSlew( dtvec, tvec, wvform );
% [ wvform_acc, dtvec_acc, tvec_acc ] = calculateSlew( dtvec_slew, tvec_slew, wvform_slew );

end