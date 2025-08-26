function [ Aspatial, bspatial ] = generateSPINSSTAMatrices( tvec, dtvec, s, staSt, opt )

%% Get k-space locations
% Get initial k-space locations
k_noadj = kxyz_fn( tvec - staSt.initSlewTime, s );

% Adjust for the final slew
k_adjfinalslew = k_noadj;
k_adjfinalslew( :, staSt.spins_int_idx ) = k_adjfinalslew( :, staSt.spins_int_idx )...
    - kxyz_fn( s.T-staSt.finalSlewTime, s )...
    - ( opt.gyro * 0.5 * staSt.finalSlewTime ) * Gxyz_fn( s.T-staSt.finalSlewTime, s );
k_adjfinalslew( :, staSt.final_slew_i:staSt.final_slew_f ) = ...
    - ( ( opt.gyro * Gxyz_fn( s.T - staSt.finalSlewTime, s ) ) / ( 2 * staSt.finalSlewTime ) ) *...
    ( ( tvec( staSt.final_slew_i:staSt.final_slew_f ) - s.T ).^2 );
k_adjfinalslew = k_adjfinalslew( :, staSt.spins_int_i:staSt.final_slew_f  );

% Adjust for initial slew
G0 = Gxyz_fn( 0, s );
deltak_initslew = ( opt.gyro / ( 2 * staSt.initSlewTime ) ) * G0 * ( tvec( staSt.init_slew_i:staSt.init_slew_f ).^2 - staSt.initSlewTime^2 );
k_initslew = deltak_initslew + kxyz_fn( 0, s );

k = [ k_initslew, k_adjfinalslew ];
k_repmat = repmat( k, [ 1, opt.numXYCoils ] );

tvec_repmat = repmat( tvec, [ 1, opt.numXYCoils ] );
dtvec_repmat = repmat( dtvec, [ 1, opt.numXYCoils ] );

Aspatial = ...
    ( 1j * ( opt.gyro * ( dtvec_repmat ) ) .*...
        ( exp( 1j * ( ( opt.pos * k_repmat ) + opt.db0 * ( opt.gyro * ( tvec_repmat - staSt.pulseLength ) ) ) ) ) );

b1p = complex( opt.b1preal, opt.b1pimag );

for nn=1:opt.numXYCoils
    timeidxs = (nn-1) * staSt.numTimePoints + ( 1:staSt.numTimePoints );
    Aspatial( :, timeidxs ) = b1p( :, nn ) .* Aspatial( :, timeidxs );
end

if nargout > 1
    bspatial = generateComplexFlipAngleTarget( opt );
end
end