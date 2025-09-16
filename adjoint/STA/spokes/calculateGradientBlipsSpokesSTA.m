function G = calculateGradientBlipsSpokesSTA( K, opt, spokes )

numBlip = size( K, 2 ) - 1;
G = zeros( 2, numBlip );

for nn = 1:numBlip

    G( :, nn ) = (...
        ( 2 * pi ) / ( opt.gyro * ( spokes.Grad_ss_blips_const(nn) + spokes.Grad_ss_blips_slew(nn) ) ) )...
        * ( K( :, (nn+1) ) - K( :, nn ) ) ;
end

G = [ G; transpose( spokes.Grad_ss_blipVals ) ];

end