function ky = ky_fn( t, s )
ky = kr_fn( t, s ) .* sin( ktheta_fn( t, s ) ) .* sin( kphi_fn( t, s ) );
end