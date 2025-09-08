function [ AE, bE ] = catOMPMatrices( AE, bE, kTP, kTTime, pulse, opt )

numK = size( K, 2 );


AEi = makeOPMColumns( complex( opt.b1preal, opt.b1pimag ), kTp, kTTime, opt, pulse );

AE = [ AEi, AE ];

bphase = complex( opt.Mtarg( :, 1 ), opt.Mtarg( :, 2 ) ) ./ sqrt( opt.Mtarg( :, 1 ).^2 + opt.Mtarg( :, 2 ).^2 );
bmag = acos( opt.Mtarg( :, 3 ) ./ opt.M0( :, 3 ) );
bE = bmag .* bphase;

end