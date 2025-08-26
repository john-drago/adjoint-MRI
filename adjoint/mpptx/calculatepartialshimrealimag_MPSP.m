function [ partamagpartareal, partaphpartareal, partamagpartaimag, partaphpartaimag ] =...
    calculatepartialshimrealimag_MPSP( wv )

partamagpartareal =  wv.shimreal_MPSP ./ wv.shim_mag_MPSP;
partaphpartareal = - wv.shimimag_MPSP ./ ( wv.shim_mag_MPSP ).^2 ;
partamagpartaimag =  wv.shimimag_MPSP ./ wv.shim_mag_MPSP;
partaphpartaimag =   wv.shimreal_MPSP ./ ( wv.shim_mag_MPSP ).^2 ;

zeromagidx = wv.shim_mag_MPSP < eps( 1e2 );

if any( zeromagidx )
    zeromagnum = sum( zeromagidx );

    partamagpartareal( zeromagidx ) = ones( zeromagnum, 1 );
    partaphpartareal( zeromagidx ) = zeros( zeromagnum, 1 );
    partamagpartaimag( zeromagidx ) = zeros( zeromagnum, 1 );
    partaphpartaimag( zeromagidx ) = ones( zeromagnum, 1 );
end
end
