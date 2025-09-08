function kxyz = kxyz_fn( t, s )
t = t( : ).';

kr = kr_fn( t, s );
ktheta = ktheta_fn( t, s );
kphi = kphi_fn( t, s );

kx = kr .* sin( ktheta ) .* cos( kphi );
ky = kr .* sin( ktheta ) .* sin( kphi );
kz = kr .* cos( ktheta );

kxyz = [...
    kx;...
    ky;...
    kz;...
    ];

end