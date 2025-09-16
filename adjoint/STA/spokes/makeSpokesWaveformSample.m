function wv = makeSpokesWaveformSample( RF0init, G0init, spokes )

%% Initialize arrays
tvec = [];
dwxyvec = [];
brealphasor = [];
bimagphasor = [];
Grad = [];

%% Prepare
% dt_tol = spokes.dt/10;

%% Loop across the spokes
for nn = 1:spokes.numSpokes
    if ( nn > 1 ) && ( nn <= ( spokes.numSpokes ) ) % handle the Grad blip before wave

        tdom_blip = spokes.tdom_blips( nn-1, : );
        tvec_blip = spokes.tvec_blips{ nn-1 };
        dwxyvec_blip = spokes.dwxyvec_blips{ nn-1 };
        
        Grad_blip = repmat( G0init( :, nn-1 ), [ 1, length( tvec_blip ) ] );
        Grad_blip( :, 1:spokes.Grad_ss_blips_slew_i( nn-1, 1 ) ) = ( G0init( :, nn-1 ) / spokes.Grad_ss_blips_slew( nn-1 ) )...
                * ( tvec_blip( 1:spokes.Grad_ss_blips_slew_i( nn-1, 1 ) ) - tdom_blip( 1 ) );
        Grad_blip( :, spokes.Grad_ss_blips_slew_f( nn-1, 1 ):end ) = -( G0init( :, nn-1 ) / spokes.Grad_ss_blips_slew( nn-1 ) )...
                * ( tvec_blip( spokes.Grad_ss_blips_slew_f( nn-1, 1 ):end ) - tdom_blip( 2 ) );

        Grad = [ Grad, Grad_blip ]; %#ok

        tvec = [ tvec, tvec_blip ]; %#ok

        dwxyvec = [ dwxyvec, dwxyvec_blip ]; %#ok

        brealphasor = [ brealphasor, zeros( size( brealphasor, 1 ), length( tvec_blip ) ) ]; %#ok
        bimagphasor = [ bimagphasor, zeros( size( bimagphasor, 1 ), length( tvec_blip ) ) ]; %#ok
        
    end
    
    tvec = [ tvec, spokes.tdom_spokes( nn, 1 ), spokes.tvec_spokes{ nn }, spokes.tdom_spokes( nn, 2 ) ]; %#ok

    dwxyvec = [ dwxyvec, 2*pi * spokes.RF_freq( nn ), spokes.dwxyvec_spokes{ nn }, 2*pi * spokes.RF_freq( nn ) ]; %#ok

    RF_spoke = RF0init( :, nn ) * spokes.RF_spokes{ nn };

    brealphasor = [ brealphasor, zeros( size( RF_spoke, 1 ), 1 ), real( RF_spoke ), zeros( size( RF_spoke, 1 ), 1 ) ]; %#ok
    bimagphasor = [ bimagphasor, zeros( size( RF_spoke, 1 ), 1 ), imag( RF_spoke ), zeros( size( RF_spoke, 1 ), 1 ) ]; %#ok
        
    Grad_spoke = zeros( 3, length( spokes.tvec_spokes{ nn } ) );
    Grad_spoke( 3, : ) = spokes.Grad_ss_spokes{ nn };

    Grad = [ Grad, zeros( 3, 1 ), Grad_spoke, zeros( 3, 1 ) ]; %#ok

    if nn == spokes.numSpokes % add in the rephaser

        tvec = [ tvec, spokes.tvec_rephase, spokes.pulseLength ]; %#ok

        dwxyvec = [ dwxyvec, spokes.dwxyvec_rephase, 2*pi * spokes.RF_freq( end ) ]; %#ok

        brealphasor = [ brealphasor, zeros( size( brealphasor, 1 ), length( spokes.tvec_rephase )+1 ) ]; %#ok
        bimagphasor = [ bimagphasor, zeros( size( bimagphasor, 1 ), length( spokes.tvec_rephase )+1 ) ]; %#ok

        Grad_rephase = zeros( 3, length( spokes.Grad_ss_rephase ) );
        Grad_rephase( 3, : ) = spokes.Grad_ss_rephase;

        Grad = [ Grad, Grad_rephase, zeros( 3, 1 ) ]; %#ok
    end

end

%% Add to wv struct
wv.tvec = tvec;
wv.dwxyvec = dwxyvec;
wv.brealphasor = brealphasor;
wv.bimagphasor = bimagphasor;
wv.breal = brealphasor;
wv.bimag = bimagphasor;
wv.RFphasor = complex( brealphasor, bimagphasor );
wv.RF = wv.RFphasor;
wv.Grad = Grad;
wv.Shim = zeros( size( Grad( 1, : ) ) );
wv.pulseLength = spokes.pulseLength;

end