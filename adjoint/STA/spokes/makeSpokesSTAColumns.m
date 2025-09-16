function spokesSTAcols = makeSpokesSTAColumns( b1p, spokesLocXY, spokeIdx, opt, spokes )

gyro = opt.gyro;

% Get the exponential contribution
k_ss = spokes.K_ss_spokes{ spokeIdx };
k_x = spokesLocXY( 1 ) * ones( size( k_ss ) );
k_y = spokesLocXY( 2 ) * ones( size( k_ss ) );

kArray = [ k_x; k_y; k_ss ];

db0 = opt.db0;
% db0 = db0 + ( -1 * 2*pi / opt.gyro ) * spokes.RF_freq( spokeIdx ); % add here if we want to change rotating frame

expArg = 2*pi * ( opt.pos * kArray ) + db0 * ( gyro * ( spokes.tvec_spokes{ spokeIdx } - spokes.pulseLength ) );

% Get the RF contribution
% RFcont = ( ( spokes.RF_spokes{ spokeIdx } ) .* spokes.dtvec_spokes{ spokeIdx } ); % keep here if we want to change rotating frame
RFcont = ( ( spokes.RF_spokes{ spokeIdx } ) .* spokes.dtvec_spokes{ spokeIdx } )...
    .* exp( 1j * ( -2*pi * spokes.RF_freq( spokeIdx ) * spokes.tvec_spokes{ spokeIdx } ) ); % how we would simulate everything in the Larmor frame

% put together
STAint = sum( RFcont .* exp( 1j * expArg ), 2 );

% STAint = sum(...
%     ( ( spokes.RF_spokes{ spokeIdx } ) .* spokes.dtvec_spokes{ spokeIdx } ) .*...
%     exp( 1j * 2*pi * ( opt.pos * kArray ) ) .*...
%     exp( 1j * opt.db0 * ( gyro * ( spokes.tvec_spokes{ spokeIdx } - spokes.pulseLength ) ) ) .*...
%     exp( 1j * ( -2*pi*spokes.modFreq * ( spokes.tvec_spokes{ spokeIdx } - spokes.pulseLength ) ) ), 2 );

% STAint = sum(...
%     ( ( spokes.RF_spokes{ spokeIdx } .* spokes.RF_mod{ spokeIdx } ) .* spokes.dtvec_spokes{ spokeIdx } ) .*...
%     exp( 1j * expArg ), 2 );

spokesSTAcols = ( opt.gyro * 1j ) * ( ( b1p ) .* STAint );

end