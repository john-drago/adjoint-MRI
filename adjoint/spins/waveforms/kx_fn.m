function kx = kx_fn( t, s )
kx = kr_fn( t, s ) .* sin( ktheta_fn( t, s ) ) .* cos( kphi_fn( t, s ) );
end