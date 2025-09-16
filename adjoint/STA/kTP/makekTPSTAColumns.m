function kTPSTAcols = makekTPSTAColumns( b1p, kTp, kTTime, opt, pulse )

gyro = opt.gyro;
pulseLength = pulse.length;
RFLength = pulse.RFLength;
minRFSlewTime = pulse.minRFSlewTime;

kTPSTAcols = ( b1p )  .* ...
        ( 1j * ( gyro * ( RFLength - minRFSlewTime ) ) *...
        ( exp( 1j * ( 2*pi * ( opt.pos * kTp ) + opt.db0 * ( gyro * ( kTTime - pulseLength ) ) ) ) ) );

end