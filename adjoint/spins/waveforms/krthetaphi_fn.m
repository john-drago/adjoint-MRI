function krthetaphi = krthetaphi_fn( t, s )
t = t( : ).';

kr = kr_fn( t, s );
ktheta = ktheta_fn( t, s );
kphi = kphi_fn( t, s );

krthetaphi = [...
    kr;...
    ktheta;...
    kphi;...
    ];

end