function G = calculateGradientBlipskTPSTA( K, opt, pulse )

numBlip = size( K, 2 ) - 1;
G = zeros( 3, numBlip );

for nn = 1:numBlip
    G( :, nn ) = ( ( 2 * 2 * pi ) / ( opt.gyro * pulse.blipLength ) ) * ( K( :, (nn+1) ) - K( :, nn ) ) ;
end

end