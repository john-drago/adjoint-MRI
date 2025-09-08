function [ partialfpartialbrealfourier, partialfpartialbimagfourier ] = calcpartialBpartialRFFourier(...
    partialfpartialBx, partialfpartialBy, opt, wv, nn  )

FBRF_nn = wv.FBRF( nn, : );

dbxdbreal = +real( FBRF_nn );
dbydbreal = +imag( FBRF_nn );
dbxdbimag = -imag( FBRF_nn );
dbydbimag = +real( FBRF_nn );

b1preal_rshp = reshape( +wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] );
b1pimag_rshp = reshape( +wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] );

dBxdbreal = reshape( +b1preal_rshp .* dbxdbreal + -b1pimag_rshp .* dbydbreal, [ opt.numPos, 1, opt.numFourier_RF*opt.numXYCoils ] );
dBydbreal = reshape( +b1pimag_rshp .* dbxdbreal + +b1preal_rshp .* dbydbreal, [ opt.numPos, 1, opt.numFourier_RF*opt.numXYCoils ] );
dBxdbimag = reshape( +b1preal_rshp .* dbxdbimag + -b1pimag_rshp .* dbydbimag, [ opt.numPos, 1, opt.numFourier_RF*opt.numXYCoils ] );
dBydbimag = reshape( +b1pimag_rshp .* dbxdbimag + +b1preal_rshp .* dbydbimag, [ opt.numPos, 1, opt.numFourier_RF*opt.numXYCoils ] );

partialfpartialbrealfourier = ...
    +partialfpartialBx .* dBxdbreal +...
    +partialfpartialBy .* dBydbreal;

partialfpartialbimagfourier =...
    +partialfpartialBx .* dBxdbimag +...
    +partialfpartialBy .* dBydbimag;

end