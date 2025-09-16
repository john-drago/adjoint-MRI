function [ AE, bE ] = generatekTPSTAMatrices( K, pulse, opt )

numK = size( K, 2 );

AE = zeros( size( opt.Mtarg, 1 ), opt.numXYCoils * numK );

pulseLength = pulse.length;
RFLength = pulse.RFLength;
% minRFSlewTime = pulse.minRFSlewTime;
blipLength = pulse.blipLength;

for kk = uint32( numK : -1 : 1  )
    
    kTidx = double( numK - kk + 1 );
    kTTime = pulseLength...
        - ( RFLength/2 )...
        - ( kTidx - 1 ) * ( RFLength ) - ( kTidx - 1 ) * blipLength;

    Acolidx = uint32( 1:opt.numXYCoils ) + ( kk - 1 ) * opt.numXYCoils;
    AE( :, Acolidx ) = makekTPSTAColumns(...
        complex( opt.b1preal, opt.b1pimag ), K( :, kk ), kTTime, opt, pulse );
end

if nargout > 1
    bE = generateComplexFlipAngleTarget( opt );
end
end