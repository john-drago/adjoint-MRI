function [ AE, bE ] = generateSpokesSTAMatrices( K, spokes, opt )

numK = size( K, 2 );

AE = zeros( size( opt.Mtarg, 1 ), opt.numXYCoils * numK );

for spokeIdx = uint32( spokes.numSpokes : -1 : ( spokes.numSpokes - numK + 1 ) )
    
    spokeCtr = spokes.numSpokes - spokeIdx + 1;

    kkDown = double( numK - spokeCtr + 1 );

    Acolidx = uint32( 1:opt.numXYCoils ) + ( kkDown - 1 ) * opt.numXYCoils;
    AE( :, Acolidx ) = makeSpokesSTAColumns(...
        complex( opt.b1preal, opt.b1pimag ), K( :, kkDown ), spokeIdx, opt, spokes );
end

if nargout > 1
    bE = generateComplexFlipAngleTarget( opt );
end
end