function bspatial = generateComplexFlipAngleTarget( opt )

% Flip angle (should work better for large flip angles)
% see Boult and Hoult, MRM, 2012. doi: 10.1002/mrm.23270

Mtargabsxy = sqrt( opt.Mtarg( :, 1 ).^2 + opt.Mtarg( :, 2 ).^2 );
bphase = complex( opt.Mtarg( :, 1 ), opt.Mtarg( :, 2 ) ) ./ Mtargabsxy;
bphase( Mtargabsxy == 0 ) = 0;
bmag = acos( opt.Mtarg( :, 3 ) ./ opt.M0( :, 3 ) );
bspatial = bmag .* bphase;

% Mtarg
% bspatial = abs( complex( opt.Mtarg( :, 1 ), opt.Mtarg( :, 2 ) ) );
end