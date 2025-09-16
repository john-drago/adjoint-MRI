function [ wvform_slew, dtvec_slew, tvec_slew ] = calculateSlew( dtvec, tvec, wvform )
% This function will calculate the slew rate of a set of waveforms, where
% each row is a distinct channel's waveform.

% expand wvform, tvec and dtvec
dtvec = transpose( dtvec( : ) );
tvec = transpose( tvec( : ) );

if size( wvform, 2 ) == length( tvec )
elseif size( wvform, 1 ) == length( tvec )
    wvform = transpose( wvform );
else
    error( "Size of waveform not compatible for interpolation." )
end

wvform_exp = [ zeros( size( wvform, 1 ), 1 ), wvform, zeros( size( wvform, 1 ), 1 ) ];
tvec_exp = [ 0, tvec, tvec(end)+dtvec(end)/2 ];

% calculate forward difference
dim = 2;
n = 1;

dtvec_slew = diff( tvec_exp, n, dim );
tvec_slew = [ 0, tvec(1:(end))] + dtvec_slew / 2;

dtvec_slew_div = dtvec_slew;
dtvec_slew_div( 1 ) = dtvec(1);
dtvec_slew_div( end ) = dtvec(end);
wvform_slew = diff( wvform_exp, n, dim ) ./ dtvec_slew_div;
% tvec_slew = [ 0, cumsum( dtvec_slew(1:(end-1)) )] + dtvec_slew / 2;


% assign zeros where we divide by 0
wvform_slew( isnan( wvform_slew ) ) = 0;
end