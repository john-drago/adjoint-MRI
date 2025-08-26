function dkxyzdparam = dkxyzdparam_fn( t, s )

%% Initialize derivatives
t = t( : ).';

indenom = indenom_fn( t, s );
denom = denom_fn( t, s );
kr = kr_fn( t, s );
ktheta = ktheta_fn( t, s );
kphi = kphi_fn( t, s );
sinktheta = sin( ktheta );
cosktheta = cos( ktheta );
sinkphi = sin( kphi );
coskphi = cos( kphi );
% dkrdt = dkrdt_fn( t, s );
% dkthetadt = dkthetadt_fn( t, s );
% dkphidt = dkphidt_fn( t, s );
% dsinkthetadt = dkthetadt .* cos( ktheta );
% dcoskthetadt = -dkthetadt .* sin( ktheta );
% dsinkphidt = dkphidt .* coskphi;
% dcoskphidt = -dkphidt .* sinkphi;

% kmax
dkrdkmax = 1 ./ denom;

% a
dkrda = ( -s.kmax * ( ( t/s.T - s.b ) .* indenom ) ) ./ ( denom.^2 );

% b
dkrdb = ( ( s.kmax * s.a ) .* indenom ) ./ ( denom.^2 );

% u
dsinkthetadu = t .* cos( ktheta );
dcoskthetadu = -t .* sin( ktheta );

% v
dsinkphidv = t .* coskphi;
dcoskphidv = -t .* sinkphi;

dkxyzdparam = zeros( 3, length( t ), 5 );

%% partial kmax

dkxdkmax = dkrdkmax .* sinktheta .* coskphi;
dkydkmax = dkrdkmax .* sinktheta .* sinkphi;
dkzdkmax = dkrdkmax .* cosktheta;

dkxyzdkmax = [...
    dkxdkmax;...
    dkydkmax;...
    dkzdkmax;...
    ];

dkxyzdparam( :, :, 1 ) = dkxyzdkmax;

%% partial a

dkxda = dkrda .* sinktheta .* coskphi;
dkyda = dkrda .* sinktheta .* sinkphi;
dkzda = dkrda .* cosktheta;

dkxyzda = [...
    dkxda;...
    dkyda;...
    dkzda;...
    ];

dkxyzdparam( :, :, 2 ) = dkxyzda;

%% partial b

dkxdb = dkrdb .* sinktheta .* coskphi;
dkydb = dkrdb .* sinktheta .* sinkphi;
dkzdb = dkrdb .* cosktheta;

dkxyzdb = [...
    dkxdb;...
    dkydb;...
    dkzdb;...
    ];

dkxyzdparam( :, :, 3 ) = dkxyzdb;

%% partial u

dkxdu = kr .* dsinkthetadu .* coskphi;
dkydu = kr .* dsinkthetadu .* sinkphi;
dkzdu = kr .* dcoskthetadu;

dkxyzdu = [...
    dkxdu;...
    dkydu;...
    dkzdu;...
    ];

dkxyzdparam( :, :, 4 ) = dkxyzdu;

%% partial v

dkxdv = kr .* sinktheta .* dcoskphidv;
dkydv = kr .* sinktheta .* dsinkphidv;
dkzdv = zeros( size( t ) );

dkxyzdv = [...
    dkxdv;...
    dkydv;...
    dkzdv;...
    ];

dkxyzdparam( :, :, 5 ) = dkxyzdv;

end