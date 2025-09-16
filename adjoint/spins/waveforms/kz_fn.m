function kz = kz_fn( t, s )
kz = kr_fn( t, s ) .* cos( ktheta_fn( t, s ) );
end