function spokes = determineSpokesTimings( spokesTim )

numSpokes = spokesTim.numSpokes;
sliceDirection = spokesTim.sliceDirection;
sliceBounds = spokesTim.sliceBounds;
sliceLocation = spokesTim.sliceLocation;

%% Determine pulse timings
if ~isfield( spokesTim, 'dt' )
    dt = 2.5e-6;
else
    dt = spokesTim.dt;
end

dttol = dt/10;
numDigAfterDecimal = 10;

if numSpokes > 1
    pulseLength = spokesTim.rephaseLength + spokesTim.centralSpokeLength + spokesTim.lastBlipLength + ...
        ( spokesTim.numSpokes - 2 ) * spokesTim.spokeBlipLength + ( spokesTim.numSpokes - 1 ) * ( spokesTim.nonCentralSpokeLength );
else
    pulseLength = spokesTim.rephaseLength + spokesTim.centralSpokeLength;
end

% Start backwards and assign timing and gradient values

% rephase period
% assume negative gradient is played here
tdom_rephase = [ pulseLength-spokesTim.rephaseLength, pulseLength ];
tvec_rephase = (tdom_rephase( 1 ) + dt/2): dt : tdom_rephase(end);
dtvec_rephase = dt * ones( size( tvec_rephase ) );
dwxyvec_rephase = -2 * pi * spokesTim.modFreq * ones( size( tvec_rephase ) );

Grad_ss_rephase = zeros( 1, length( tvec_rephase ) );
Grad_ss_rephase_slew_i = find( tvec_rephase <= ( tdom_rephase( 1 ) + spokesTim.rephaseSlewTime + dttol ), 1, 'last' );
Grad_ss_rephase_slew_f = find( tvec_rephase >= ( tdom_rephase( end ) - spokesTim.rephaseSlewTime - dttol  ), 1, 'first' );
Grad_ss_rephase_const_idx = [ Grad_ss_rephase_slew_i+1, Grad_ss_rephase_slew_f-1 ];
if diff( Grad_ss_rephase_const_idx ) < 0
    Grad_ss_rephase_const_idx = zeros( 0, 2 );
end
Grad_ss_rephase( 1, : ) = -spokesTim.rephaseGradMag;
Grad_ss_rephase( 1, 1:Grad_ss_rephase_slew_i ) = ( -spokesTim.rephaseGradMag / spokesTim.rephaseSlewTime )...
    * ( tvec_rephase( 1:Grad_ss_rephase_slew_i ) - tdom_rephase( 1 ) );
Grad_ss_rephase( 1, Grad_ss_rephase_slew_f:end) = ( spokesTim.rephaseGradMag / spokesTim.rephaseSlewTime )...
    * ( tvec_rephase( Grad_ss_rephase_slew_f:end ) - tdom_rephase( end ) );
GradTimeProduct_rephase = Grad_ss_rephase( 1, : ) .* dtvec_rephase;

K_ss_rephase = zeros( size( tvec_rephase ) );
K_ss_rephase( 1, : ) = - spokesTim.gyro/( 2*pi ) * ( cumsum( GradTimeProduct_rephase, 2, 'reverse' ) +...
    - 0.5 * GradTimeProduct_rephase);
K_ss_st_previous = K_ss_rephase( 1, 1 ) - spokesTim.gyro/( 2*pi ) * 0.5 * GradTimeProduct_rephase( 1 );

% determine gradients in the blip before the central spoke
% if numSpokes > 1
%     K_ss_lastblip = - spokesTim.gyro / ( 2*pi ) * spokesTim.lastBlipGradMag * ( spokesTim.lastBlipConstantTime + spokesTim.lastBlipSlewTime );
% end

% spoke periods
tdom_spokes = zeros( numSpokes, 2 );
tvec_spokes = cell( numSpokes, 1 );
dtvec_spokes = cell( numSpokes, 1 );
dwxyvec_spokes = cell( numSpokes, 1 );
length_spokes = zeros( numSpokes, 1 );
Grad_ss_spokes = cell( numSpokes, 1 ); 
Grad_ss_mag_spokes = zeros( numSpokes, 1 );
Grad_ss_spokes_slew = zeros( numSpokes, 1 );
Grad_ss_spokes_const = zeros( numSpokes, 1 );
Grad_ss_spokes_length = zeros( numSpokes, 1 );
Grad_ss_spokes_slew_i = zeros( numSpokes, 1 );
Grad_ss_spokes_slew_f = zeros( numSpokes, 1 );
Grad_ss_spokes_const_idx = cell( numSpokes, 1 );
K_ss_spokes = cell( numSpokes, 1 );
K_ss_st_spokes = zeros( numSpokes, 1 );
K_ss_minmax_spokes = zeros( numSpokes, 2 );

Grad_ss_blipVals = zeros( numSpokes-1, 1 );
Grad_ss_blips_slew = zeros( numSpokes-1, 1 );
Grad_ss_blips_const = zeros( numSpokes-1, 1 );
Grad_ss_blips_length = zeros( numSpokes-1, 1 );
Grad_ss_blips = cell( numSpokes-1, 1 ); 
Grad_ss_blips_slew_i = zeros( numSpokes-1, 1 );
Grad_ss_blips_slew_f = zeros( numSpokes-1, 1 );
Grad_ss_blips_const_idx = cell( numSpokes-1, 1 );
K_ss_blips = cell( numSpokes-1, 1 );
tdom_blips = zeros( numSpokes-1, 2 );
tvec_blips = cell( numSpokes-1, 1 );
dtvec_blips = cell( numSpokes-1, 1 );
dwxyvec_blips = cell( numSpokes-1, 1 );
length_blips = zeros( numSpokes-1, 1 );

Grad_ss = Grad_ss_rephase;
tvec_ss = tvec_rephase;
dtvec_ss = dtvec_rephase;
dwxyvec_ss = dwxyvec_rephase;
K_ss = K_ss_rephase;
RF_ss = zeros( size( tvec_rephase  ) );

RF_spokes = cell( numSpokes, 1 );
RF_freq = zeros( numSpokes, 1 );
RF_BW = zeros( numSpokes, 1 );
TBW_spokes = zeros( numSpokes, 1 );

for nn = (numSpokes):-1:1
    spokeIdx = numSpokes - nn + 1;

    if mod( spokeIdx, 2 ) == 1
        gradSign = +1;
    else
        gradSign = -1;
    end

    if nn == numSpokes
        tdom_spokes( nn, : ) = round( ...
            [...
            pulseLength - spokesTim.rephaseLength - spokesTim.centralSpokeLength,...
            pulseLength - spokesTim.rephaseLength ], numDigAfterDecimal );
        Grad_ss_spokes_slew( nn ) = spokesTim.centralSpokeSlewTime;
        Grad_ss_spokes_const( nn ) = spokesTim.centralSpokeConstantTime;
        Grad_ss_spokes_length( nn ) = spokesTim.centralSpokeLength;
        Grad_ss_mag_spokes( nn ) = gradSign * spokesTim.centralSpokeGradMag;
        TBW_spokes( nn ) = spokesTim.centralSpokeTBW;
        length_spokes( nn ) = spokesTim.centralSpokeLength;

    elseif nn == (numSpokes - 1)

        tdom_spokes( nn, : ) = round(...
            [...
            tdom_spokes( nn+1, 1 ) - spokesTim.nonCentralSpokeLength - spokesTim.lastBlipLength,...
            tdom_spokes( nn+1, 1 ) - spokesTim.lastBlipLength ], numDigAfterDecimal );
        
        % K_ss_st_previous = K_ss_st_previous + K_ss_lastblip;

        Grad_ss_spokes_slew( nn ) = spokesTim.nonCentralSpokeSlewTime;
        Grad_ss_spokes_const( nn ) = spokesTim.nonCentralSpokeConstantTime;
        Grad_ss_spokes_length( nn ) = spokesTim.nonCentralSpokeLength;
        Grad_ss_mag_spokes( nn ) = gradSign * spokesTim.nonCentralSpokeGradMag;
        TBW_spokes( nn ) = spokesTim.nonCentralSpokeTBW;
        length_spokes( nn ) = spokesTim.nonCentralSpokeLength;

        Grad_ss_blipVals( nn ) = spokesTim.lastBlipGradMag;
        Grad_ss_blips_slew( nn ) = spokesTim.lastBlipSlewTime;
        Grad_ss_blips_const( nn ) = spokesTim.lastBlipConstantTime;
        Grad_ss_blips_length( nn ) = spokesTim.lastBlipLength;

        length_blips( nn ) = spokesTim.lastBlipLength;
        tdom_blips( nn, : ) = round( [...
            tdom_spokes( nn+1, 1 ) - spokesTim.lastBlipLength,...
            tdom_spokes( nn+1, 1 ) ], numDigAfterDecimal );

    else

        tdom_spokes( nn, : ) = round( ...
            [...
            tdom_spokes( nn+1, 1 ) - spokesTim.nonCentralSpokeLength - spokesTim.spokeBlipLength,...
            tdom_spokes( nn+1, 1 ) - spokesTim.spokeBlipLength ], numDigAfterDecimal );
        Grad_ss_spokes_slew( nn ) = spokesTim.nonCentralSpokeSlewTime;
        Grad_ss_spokes_const( nn ) = spokesTim.nonCentralSpokeConstantTime;
        Grad_ss_spokes_length( nn ) = spokesTim.nonCentralSpokeLength;
        Grad_ss_mag_spokes( nn ) = gradSign * spokesTim.nonCentralSpokeGradMag;
        TBW_spokes( nn ) = spokesTim.nonCentralSpokeTBW;
        length_spokes( nn ) = spokesTim.nonCentralSpokeLength;

        Grad_ss_blips_slew( nn ) = spokesTim.spokeBlipLength / 2;
        Grad_ss_blips_const( nn ) = 0;
        Grad_ss_blips_length( nn ) = spokesTim.spokeBlipLength;
        length_blips( nn ) = spokesTim.spokeBlipLength;
        tdom_blips( nn, : ) = round( [...
            tdom_spokes( nn+1, 1 ) - spokesTim.spokeBlipLength,...
            tdom_spokes( nn+1, 1 ) ], numDigAfterDecimal );
    end

    if nn < numSpokes
        tvec_blips{ nn } = round( ( tdom_blips( nn, 1 ) + dt/2 ) : dt : ( tdom_blips( nn, 2 ) ), numDigAfterDecimal );
        dtvec_blips{ nn } = round( dt * ones( size( tvec_blips{ nn } ) ), numDigAfterDecimal );
        dwxyvec_blips{ nn } = ( -2 * pi * gradSign * spokesTim.modFreq ) * ones( size( tvec_blips{ nn } ) );

        % tcenter_blip = sum( tdom_blips( nn, : ) ) / 2;
        % Grad_ss_blips_rise_idx{ nn } = find( tvec_blips{ nn } <= ( tcenter_blip + dttol ) );
        % Grad_ss_blips_fall_idx{ nn } = find( tvec_blips{ nn } >= ( tcenter_blip - dttol ) );

        Grad_ss_blips{ nn } = zeros( size( tvec_blips{ nn } ) );
        Grad_ss_blips{ nn }( 1, : ) = Grad_ss_blipVals( nn );
        
        blips_slew_i = find(...
            tvec_blips{ nn } <= ( tdom_blips( nn, 1 ) + Grad_ss_blips_slew( nn ) + dttol ), 1, 'last' );
        blips_slew_f = find(...
            tvec_blips{ nn } >= ( tdom_blips( nn, 2 ) - Grad_ss_blips_slew( nn ) - dttol  ), 1, 'first' );

        if ~isempty( blips_slew_i )
            Grad_ss_blips_slew_i( nn, 1 ) = blips_slew_i;
            Grad_ss_blips{ nn }( 1, 1:Grad_ss_blips_slew_i( nn, 1 ) ) = ( Grad_ss_blipVals( nn ) / Grad_ss_blips_slew( nn ) )...
                * ( tvec_blips{ nn }( 1:Grad_ss_blips_slew_i( nn, 1 ) ) - tdom_blips( nn, 1 ) );
        end
        if ~isempty( blips_slew_f )
            Grad_ss_blips_slew_f( nn, 1 ) = blips_slew_f;
            Grad_ss_blips{ nn }( 1, Grad_ss_blips_slew_f( nn, 1 ):end) = -( Grad_ss_blipVals( nn ) / Grad_ss_blips_slew( nn ) )...
                * ( tvec_blips{ nn }( Grad_ss_blips_slew_f( nn, 1 ):end ) - tdom_blips( nn, 2 ) );
        end
        if ~isempty( blips_slew_i ) && ~isempty( blips_slew_f )
            Grad_ss_blips_const_idx{ nn } = [...
                Grad_ss_blips_slew_i( nn, 1 ) + 1,...
                Grad_ss_blips_slew_f( nn, 1 ) - 1 ];
        end
        if diff( Grad_ss_blips_const_idx{ nn } ) < 0
            Grad_ss_blips_const_idx{ nn } = zeros( 0, 2 );
        end

        GradTimeProduct_blips = Grad_ss_blips{ nn }( 1, : ) .* dtvec_blips{ nn };

        K_ss_blips{ nn } = zeros( size( tvec_blips{ nn } ) );
        K_ss_blips{ nn }( 1, : ) = - spokesTim.gyro/( 2*pi ) * ( cumsum( GradTimeProduct_blips, 2, 'reverse' ) +...
            - 0.5 * GradTimeProduct_blips ) + K_ss_st_previous;
        K_ss_st_blips( nn, 1 ) = K_ss_blips{ nn }( 1, 1 ) - spokesTim.gyro/( 2*pi ) * 0.5 * GradTimeProduct_blips( 1 );
        K_ss_minmax_blips( nn, 1 ) = K_ss_st_previous;
        K_ss_st_previous = K_ss_blips{ nn }( 1, 1 ) - spokesTim.gyro/( 2*pi ) * 0.5 * GradTimeProduct_blips( 1 );
        K_ss_minmax_blips( nn, 2 ) = K_ss_st_previous;

        Grad_ss = [ Grad_ss_blips{ nn }, Grad_ss ]; %#ok
        K_ss = [ K_ss_blips{ nn }, K_ss ]; %#ok
        tvec_ss = [ tvec_blips{ nn }, tvec_ss ]; %#ok
        dtvec_ss = [ dtvec_blips{ nn }, dtvec_ss ]; %#ok
        dwxyvec_ss = [ dwxyvec_blips{ nn }, dwxyvec_ss ]; %#ok
        RF_ss = [ zeros( size( tvec_blips{ nn } ) ), RF_ss ]; %#ok

    end

    tvec_spokes{ nn } = round( (tdom_spokes( nn, 1 ) + dt/2): dt : tdom_spokes( nn, 2 ), numDigAfterDecimal );
    dtvec_spokes{ nn } = round( dt * ones( size( tvec_spokes{ nn } ) ), numDigAfterDecimal );
    dwxyvec_spokes{ nn } = ( -2 * pi * gradSign * spokesTim.modFreq ) * ones( size( tvec_spokes{ nn } ) );

    Grad_ss_spokes{ nn } = zeros( 1, length( tvec_spokes{ nn } ) );
    Grad_ss_spokes_slew_i( nn, 1 ) = find(...
        tvec_spokes{ nn } <= ( tdom_spokes( nn, 1 ) + Grad_ss_spokes_slew( nn ) + dttol ), 1, 'last' );
    Grad_ss_spokes_slew_f( nn, 1 ) = find(...
        tvec_spokes{ nn } >= ( tdom_spokes( nn, 2 ) - Grad_ss_spokes_slew( nn ) - dttol  ), 1, 'first' );
    Grad_ss_spokes_const_idx{ nn } = [...
        Grad_ss_spokes_slew_i( nn, 1 ) + 1,...
        Grad_ss_spokes_slew_f( nn, 1 ) - 1 ];
    if diff( Grad_ss_spokes_const_idx{ nn } ) < 0
        Grad_ss_spokes_const_idx{ nn } = zeros( 0, 2 );
    end
    Grad_ss_spokes{ nn }( 1, : ) = Grad_ss_mag_spokes( nn );
    Grad_ss_spokes{ nn }( 1, 1:Grad_ss_spokes_slew_i( nn, 1 ) ) = ( Grad_ss_mag_spokes( nn ) / Grad_ss_spokes_slew( nn ) )...
        * ( tvec_spokes{ nn }( 1:Grad_ss_spokes_slew_i( nn, 1 ) ) - tdom_spokes( nn, 1 ) );
    Grad_ss_spokes{ nn }( 1, Grad_ss_spokes_slew_f( nn, 1 ):end) = -( Grad_ss_mag_spokes( nn ) / Grad_ss_spokes_slew( nn ) )...
        * ( tvec_spokes{ nn }( Grad_ss_spokes_slew_f( nn, 1 ):end ) - tdom_spokes( nn, 2 ) );
    GradTimeProduct_spokes = Grad_ss_spokes{ nn }( 1, : ) .* dtvec_spokes{ nn };

    K_ss_spokes{ nn } = zeros( size( tvec_spokes{ nn } ) );
    K_ss_spokes{ nn }( 1, : ) = - spokesTim.gyro/( 2*pi ) * ( cumsum( GradTimeProduct_spokes, 2, 'reverse' ) +...
        - 0.5 * GradTimeProduct_spokes ) + K_ss_st_previous;
    K_ss_st_spokes( nn, 1 ) = K_ss_spokes{ nn }( 1, 1 ) - spokesTim.gyro/( 2*pi ) * 0.5 * GradTimeProduct_spokes( 1 );
    K_ss_minmax_spokes( nn, 1 ) = K_ss_st_previous;
    K_ss_st_previous = K_ss_spokes{ nn }( 1, 1 ) - spokesTim.gyro/( 2*pi ) * 0.5 * GradTimeProduct_spokes( 1 );
    K_ss_minmax_spokes( nn, 2 ) = K_ss_st_previous;

    K_ss_minmax_spokes( nn, : ) = round( K_ss_minmax_spokes( nn, : ), numDigAfterDecimal );

    sc_K_ss = K_ss_spokes{ nn } ;

    RF_BW( nn, 1 ) = TBW_spokes( nn, 1 ) / length_spokes( nn, 1 ); % believe that this is two-sided bandwidth
    RF_spokes{ nn } = ...
        sinc( ( 0.5 * TBW_spokes( nn ) / K_ss_minmax_spokes( nn,2 ) ) * sc_K_ss ) .* hamming( length( sc_K_ss ) - 1 );
    RF_freq( nn ) = gradSign * spokesTim.modFreq;

    Grad_ss = [ Grad_ss_spokes{ nn }, Grad_ss ]; %#ok
    K_ss = [ K_ss_spokes{ nn }, K_ss ]; %#ok
    tvec_ss = [ tvec_spokes{ nn }, tvec_ss ]; %#ok
    dtvec_ss = [ dtvec_spokes{ nn }, dtvec_ss ]; %#ok
    dwxyvec_ss = [ dwxyvec_spokes{ nn }, dwxyvec_ss ]; %#ok
    RF_ss = [ RF_spokes{ nn }, RF_ss ]; %#ok

end

%% Generate spokes struct
spokes = struct;
spokes.tdom_spokes = tdom_spokes;
spokes.tvec_spokes = tvec_spokes;
spokes.dtvec_spokes = dtvec_spokes;
spokes.dwxyvec_spokes = dwxyvec_spokes;

spokes.Grad_ss_spokes = Grad_ss_spokes;
spokes.Grad_ss_mag_spokes = Grad_ss_mag_spokes;
spokes.Grad_ss_spokes_slew = Grad_ss_spokes_slew;
spokes.Grad_ss_spokes_const = Grad_ss_spokes_const;

spokes.Grad_ss_spokes_length = Grad_ss_spokes_length;
spokes.Grad_ss_spokes_slew_i = Grad_ss_spokes_slew_i;
spokes.Grad_ss_spokes_slew_f = Grad_ss_spokes_slew_f;
spokes.Grad_ss_spokes_const_idx = Grad_ss_spokes_const_idx;

spokes.K_ss_spokes = K_ss_spokes;
spokes.K_ss_st_spokes = K_ss_st_spokes;
spokes.K_ss_minmax_spokes = K_ss_minmax_spokes;

spokes.RF_spokes = RF_spokes;
spokes.RF_freq = RF_freq;
spokes.RF_BW = RF_BW;
spokes.TBW_spokes = TBW_spokes;

spokes.sliceDirection = sliceDirection;
spokes.sliceLocation = sliceLocation;
spokes.sliceBounds = sliceBounds;
spokes.sliceThickness = spokesTim.sliceThickness;
spokes.modFreq = spokesTim.modFreq;
spokes.numSpokes = numSpokes;
spokes.pulseLength = pulseLength;
spokes.dt = dt;

spokes.K_ss = K_ss;
spokes.Grad_ss = Grad_ss;
spokes.RF_ss = RF_ss;
spokes.tvec_ss = tvec_ss;
spokes.dtvec_ss = dtvec_ss;
spokes.dwxyvec_ss = dwxyvec_ss;

spokes.Grad_ss_blipVals = Grad_ss_blipVals;
spokes.Grad_ss_blips = Grad_ss_blips;
spokes.Grad_ss_blips_slew_i = Grad_ss_blips_slew_i;
spokes.Grad_ss_blips_slew_f = Grad_ss_blips_slew_f;
spokes.Grad_ss_blips_const_idx = Grad_ss_blips_const_idx;
spokes.Grad_ss_blips_slew = Grad_ss_blips_slew;
spokes.Grad_ss_blips_const = Grad_ss_blips_const;
spokes.Grad_ss_blips_length = Grad_ss_blips_length;
spokes.K_ss_blips = K_ss_blips;
spokes.tdom_blips = tdom_blips;
spokes.tvec_blips = tvec_blips;
spokes.dtvec_blips = dtvec_blips;
spokes.dwxyvec_blips = dwxyvec_blips;
spokes.length_blips = length_blips;

spokes.tdom_rephase = tdom_rephase;
spokes.tvec_rephase = tvec_rephase;
spokes.dtvec_rephase = dtvec_rephase;
spokes.dwxyvec_rephase = dwxyvec_rephase;
spokes.Grad_ss_rephase = Grad_ss_rephase;
spokes.Grad_ss_rephase_slew_i = Grad_ss_rephase_slew_i;
spokes.Grad_ss_rephase_slew_f = Grad_ss_rephase_slew_f;
spokes.Grad_ss_rephase_const_idx = Grad_ss_rephase_const_idx;
spokes.K_ss_rephase = K_ss_rephase;

spokes.spokeBlipLength = spokesTim.spokeBlipLength;

spokes.centralSpokeSlewTime = spokesTim.centralSpokeSlewTime;
spokes.centralSpokeConstantTime = spokesTim.centralSpokeConstantTime;
spokes.centralSpokeLength = spokesTim.centralSpokeLength;
spokes.centralSpokeGradMag = spokesTim.centralSpokeGradMag;

spokes.lastBlipSlewTime = spokesTim.lastBlipSlewTime;
spokes.lastBlipConstantTime = spokesTim.lastBlipConstantTime;
spokes.lastBlipLength = spokesTim.lastBlipLength;
spokes.lastBlipGradMag = spokesTim.lastBlipGradMag;

spokes.rephaseSlewTime = spokesTim.rephaseSlewTime;
spokes.rephaseConstantTime = spokesTim.rephaseConstantTime;
spokes.rephaseLength = spokesTim.rephaseLength;
spokes.rephaseGradMag = spokesTim.rephaseGradMag;

spokes.nonCentralSpokeSlewTime = spokesTim.nonCentralSpokeSlewTime;
spokes.nonCentralSpokeConstantTime = spokesTim.nonCentralSpokeConstantTime;
spokes.nonCentralSpokeLength = spokesTim.nonCentralSpokeLength;
spokes.nonCentralSpokeGradMag = spokesTim.nonCentralSpokeGradMag;
spokes.gyro = spokesTim.gyro;

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function y = sinc( x )
pix = pi*x;
y = sin( pix ) ./ ( pix );
y( pix == 0 ) = 1;
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function h = hamming( N )
n = 0:1:N;
a0 = 0.53836;
a1 = 0.46164;
h = a0 - a1 * cos( (2*pi/N) * n );
end
% ----------------------------------------------------------------------- %