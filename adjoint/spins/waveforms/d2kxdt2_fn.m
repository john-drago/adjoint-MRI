function d2kxdt2 = d2kxdt2_fn( t, s )
t = t( : ).';

kr = kr_fn( t, s );
ktheta = ktheta_fn( t, s );
kphi = kphi_fn( t, s );
sinktheta = sin( ktheta );
% cosktheta = cos( ktheta );
sinkphi = sin( kphi );
coskphi = cos( kphi );
dkrdt = dkrdt_fn( t, s );
d2krdt2 = d2krdt2_fn( t, s );
dkthetadt = dkthetadt_fn( t, s );
dkphidt = dkphidt_fn( t, s );
dsinkthetadt = dkthetadt .* cos( ktheta );
d2sinkthetadt2 = ( - dkthetadt.^2 ) * sinktheta;
% dcoskthetadt = -dkthetadt .* sin( ktheta );
% d2coskthetadt2 = ( - dkthetadt.^2 ) * cosktheta;
% dsinkphidt = dkphidt .* coskphi;
% d2sinkphidt2 = ( - dkphidt.^2 ) .* sinkphi;
dcoskphidt = -dkphidt .* sinkphi;
d2coskphidt2 = ( - dkphidt.^2 ) .* coskphi;

d2kxdt2 = ...
    d2krdt2 .* sinktheta .* coskphi + dkrdt .* dsinkthetadt .* coskphi + dkrdt .* sinktheta .* dcoskphidt +...
    dkrdt .* dsinkthetadt .* coskphi + kr .* d2sinkthetadt2 .* coskphi + kr .* dsinkthetadt .* dcoskphidt +...
    dkrdt .* sinktheta .* dcoskphidt + kr .* dsinkthetadt .* dcoskphidt + kr .* sinktheta .* d2coskphidt2;

end